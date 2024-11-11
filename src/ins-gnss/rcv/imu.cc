/*------------------------------------------------------------------------------
* imu.cc : decode imu raw measurement data
*
* version : $Revision:$ $Date:$
* history : 2017/11/06  1.0  new
*-----------------------------------------------------------------------------*/
#include <navlib.h>

/* constants/macros ----------------------------------------------------------*/
#define NUMBYTES_GI310  43                      /* numbers of bytes of gi310 imu raw data */
#define MAXDIFFTIME     10.0                    /* max time difference to reset  */
#define NUMBYTES_ZHUFENG  32                      /* numbers of bytes of gi310 imu raw data */
const unsigned char gi310_head[2]={0x55,0xAA};  /* imu message header */
const unsigned char zhufeng_head[2]={0x24,0x49};  /* imu message header */
/* get fields (little-endian) ------------------------------------------------*/
#define U1(p) (*((unsigned char *)(p)))
#define I1(p) (*((char *)(p)))
#define EPSILON             2.2204460492503131e-016                 ///< epsilon == DBL_EPSILON
#define isEqual(a,b)        (fabs(a-b)<EPSILON ? true : false)      ///< a == b
/* get fields (little-endian) ------------------------------------------------*/
static unsigned short U2(unsigned char *p) {unsigned short u; memcpy(&u,p,2); return u;}
static unsigned int   U4(unsigned char *p) {unsigned int   u; memcpy(&u,p,4); return u;}
static short          I2(unsigned char *p) {short          i; memcpy(&i,p,2); return i;}
static int            I4(unsigned char *p) {int            i; memcpy(&i,p,4); return i;}
static float          R4(unsigned char *p) {float          r; memcpy(&r,p,4); return r;}
static double         R8(unsigned char *p) {double         r; memcpy(&r,p,8); return r;}
static char           C1(unsigned char *p) {char         c; memcpy(&c,p,1); return c;}
static unsigned char  UC1(unsigned char *p) {unsigned char         c; memcpy(&c,p,1); return c;}
static bool           B1(unsigned char *p){ bool         b;memcpy(&b,p,1);return b;}

/* for odometry option--------------------------------------------------------*/
static double res=2048; /* odometry resolution */
static double d=0.73;   /* odometry wheel diameter (m) */

/* set odometry resolution----------------------------------------------------*/
extern void odores(double rese) {res=rese;}
/* set odometry wheel diameter------------------------------------------------*/
extern void odod(double de) {d=de;}

/* check header --------------------------------------------------------------*/
static int chkhead(const unsigned char *buff, const unsigned char* head)
{
    return (buff[0]==head[0])&&(buff[1]==head[1]);
}
/* checksum ------------------------------------------------------------------*/
static unsigned char chksum(const unsigned char *buff, int len) {
    int i;
    unsigned char sum=0;
    for (i=0;i<len;i++)sum+=buff[i];return sum;
}
/* decode imu time------------------------------------------------------------*/
static void decode_sow_time(raw_t *raw,double *sow,int *start)
{
    static unsigned int pps=0;

    if (pps==0) {
        pps=U4(raw->buff+6); *sow=U4(raw->buff+2);
        return;
    }
    if (pps!=U4(raw->buff+6)) {(*sow)++; *start=1;}
    pps=U4(raw->buff+6);
}
/* adjust gps seconds of week and imu time------------------------------------*/
static void adjtime(raw_t* raw,const double sowi,double *sowo,double *timu,
                    int *week,unsigned int *dcc)
{
    static unsigned int imuc=0,dc=0;
    int d=0;
    
    if (imuc==0) {
        imuc=U4(raw->buff+10); *sowo=sowi;
        return;
    }
    *dcc=dc=U4(raw->buff+10)-imuc<0?UINT_MAX+U4(raw->buff+10)-imuc:U4(raw->buff+10)-imuc;
    imuc=U4(raw->buff+10);

    /* increase week */
    if ((*sowo=(int)(1.0/FREQOCXO*dc+sowi))>=604800.0) {
        *sowo-=604800.0; (*week)++;
    }
    d=U4(raw->buff+10)-U4(raw->buff+06);
    *timu=*sowo+1.0/FREQOCXO*d;
}
/* decode odometry and convert to velocity in vehicle frame------------------*/
static int decode_odo_data(raw_t *raw,double dt)
{
    static int dc;

    if (raw->imu.time.time==0||dt==0.0) {
        raw->imu.odoc=I2(raw->buff+38); return 0;
    }
    dc=I2(raw->buff+38)-raw->imu.odoc<=SHRT_MIN?
       I2(raw->buff+38)-raw->imu.odoc+USHRT_MAX:
       I2(raw->buff+38)-raw->imu.odoc>=SHRT_MAX?
       I2(raw->buff+38)-raw->imu.odoc-USHRT_MAX:
       I2(raw->buff+38)-raw->imu.odoc;

    raw->imu.odoc=I2(raw->buff+38);
    raw->imu.odo.time=raw->imu.time;
    raw->imu.odo.dt=dt;
    raw->imu.odo.dr=dc/res*PI*d;
    return 1;
}
/* decode imu data ----------------------------------------------------------*/
static int decode_imu_data(raw_t *raw)
{
    int i,week=0;
    unsigned int dc=0;
    static int start=0;
    static double sow=0.0,timu=0.0;
    gtime_t t0;

    trace(3,"decode_imu_data:\n");

    raw->imut.n=0;

    /* decode GPS sow (s) */
    decode_sow_time(raw,&sow,&start);

    /* start decode imu time */
    if (start) {
        adjtime(raw,sow,&sow,&timu,&week,&dc);
    }
    else return 0;

    /* current and precious time difference is too large */
    if (dc*1.0/FREQOCXO>MAXDIFFTIME) return 0;

    t0=gpst2time(week,timu);
    raw->imu.pps =U4(raw->buff+06);
    raw->imu.imuc=U4(raw->buff+10);

    for (i=0;i<3;i++) {
        raw->imu.gyro[i]=R4(raw->buff+14+i*4);
        raw->imu.accl[i]=R4(raw->buff+26+i*4);
    }
    decode_odo_data(raw,dc*1.0/FREQOCXO);
    raw->imu.time=t0;

    /* add imu measurement data */
    addimudata(&raw->imut,&raw->imu);
    return timu>0.0?4:0;
}
/* decode M39 IMU GI310 data ------------------------------------------------*/
static int decode_imu_m39(raw_t *raw, unsigned char data)
{
    trace(3,"decode_imu_m39:\n");

    raw->buff[raw->nbyte++]=data;

    if (raw->nbyte<2) return 0; /* synchronize frame */

    if (!chkhead(raw->buff,gi310_head)) {raw->nbyte=0; return 0;}
    if (raw->nbyte<NUMBYTES_GI310) return 0;

    if (chksum(raw->buff+2,40)==data) { /* checksum */
        raw->nbyte=0;
        return decode_imu_data(raw);
    }
    else {
        raw->nbyte=0; return -1; /* checksum fail */
    }
}

static int decode_imu_zhufeng(raw_t *raw, unsigned char data)
{
    trace(3,"decode_imu_zhufeng:\n");

    raw->buff[raw->nbyte++]=data;

    if (!raw->m_FFSHead.flag){
    if (raw->nbyte<512) return 0; /* synchronize frame */
    if (chkhead(raw->buff,zhufeng_head)) {
        // Read IMR File Head Part
        memcpy(raw->m_FFSHead.szHeader,raw->buff+0,8);
        raw->m_FFSHead.bIsIntelOrMotorola=C1(raw->buff+8);
        raw->m_FFSHead.dVersionNumber=R8(raw->buff+9);
        raw->m_FFSHead.bDeltaTheta=I4(raw->buff+17);
        raw->m_FFSHead.bDeltaVelocity=I4(raw->buff+21);
        raw->m_FFSHead.dDataRateHz=R8(raw->buff+25);
        raw->m_FFSHead.dGyrosScaleFactor=R8(raw->buff+33);
        raw->m_FFSHead.dAccelScaleFactor=R8(raw->buff+41);
        raw->m_FFSHead.iUtcOrGpsTime=I4(raw->buff+49);
        raw->m_FFSHead.iRcvTimeOrCorrTime=I4(raw->buff+53);
        raw->m_FFSHead.dTimeTagBias=R8(raw->buff+57);
        memcpy(raw->m_FFSHead.szImuName,raw->buff+65,32);
        raw->m_FFSHead.bDirValid=B1(raw->buff+97);
        raw->m_FFSHead.ucX=UC1(raw->buff+98);
        raw->m_FFSHead.ucY=UC1(raw->buff+99);
        raw->m_FFSHead.ucZ=UC1(raw->buff+100);
        memcpy(raw->m_FFSHead.szProgramName,raw->buff+101,32);
        raw->m_FFSHead.tCreate.year=I2(raw->buff+133);
        raw->m_FFSHead.tCreate.month=I2(raw->buff+135);
        raw->m_FFSHead.tCreate.day=I2(raw->buff+137);
        raw->m_FFSHead.tCreate.hour=I2(raw->buff+139);
        raw->m_FFSHead.tCreate.minute=I2(raw->buff+141);
        raw->m_FFSHead.tCreate.second=I2(raw->buff+143);
        if (raw->m_FFSHead.tCreate.year < 50)
            raw->m_FFSHead.tCreate.year += 2000;
        else
            raw->m_FFSHead.tCreate.year += 1900;
        raw->m_FFSHead.bLeverArmValid=B1(raw->buff+145);
        raw->m_FFSHead.lXoffset=I4(raw->buff+146);
        raw->m_FFSHead.lYoffset=I4(raw->buff+150);
        raw->m_FFSHead.lZoffset=I4(raw->buff+154);
        memcpy(raw->m_FFSHead.Reserved,raw->buff+158,354);
        raw->m_FFSHead.flag=1;
        raw->nbyte=0;
        return 0;
    }
    }
    if (raw->nbyte<NUMBYTES_ZHUFENG) return 0;
    raw->nbyte=0;
    int i,week=0;
    unsigned int dc=0;
    static int start=0;
    static double sow=0.0;
    gtime_t t0;
    time2gpst(raw->time,&week);
    trace(3,"decode_imu_data:\n");

    raw->imut.n=0;
    IE84_INS_type imrData;
    /* decode GPS sow (s) */
    sow=R8(raw->buff+0);
    imrData.gx= I4(raw->buff+8);
    imrData.gy= I4(raw->buff+12);
    imrData.gz= I4(raw->buff+16);
    imrData.ax= I4(raw->buff+20);
    imrData.ay= I4(raw->buff+24);
    imrData.az= I4(raw->buff+28);
    if (sow > 604800) sow -= 604800;
    // The unit of outputs (gx, gy, gz) for imr files is deg/s
    raw->insData.gt = sow;
    raw->insData.gx = imrData.gx * raw->m_FFSHead.dGyrosScaleFactor * raw->m_FFSHead.dDataRateHz;
    raw->insData.gy = imrData.gy * raw->m_FFSHead.dGyrosScaleFactor * raw->m_FFSHead.dDataRateHz;
    raw->insData.gz = imrData.gz * raw->m_FFSHead.dGyrosScaleFactor * raw->m_FFSHead.dDataRateHz;
    raw->insData.ax = imrData.ax * raw->m_FFSHead.dAccelScaleFactor * raw->m_FFSHead.dDataRateHz;
    raw->insData.ay = imrData.ay * raw->m_FFSHead.dAccelScaleFactor * raw->m_FFSHead.dDataRateHz;
    raw->insData.az = imrData.az * raw->m_FFSHead.dAccelScaleFactor * raw->m_FFSHead.dDataRateHz;
    // Correct the time offset
    raw->insData.gt -= (raw->m_FFSHead.dTimeTagBias * 0.001);

    if (raw->insData.gt > raw->insData.gtPre)
    {  raw->imu.time= gpst2time(week,raw->insData.gt);
       raw->imu.gyro[0] = raw->insData.gx;raw->imu.gyro[1] = raw->insData.gy;raw->imu.gyro[2] = raw->insData.gz;
       raw->imu.accl[0] = raw->insData.ax;raw->imu.accl[1] = raw->insData.ay;raw->imu.accl[2] = raw->insData.az;
       raw->insData.gtPre = raw->insData.gt;
    }
    else if (isEqual(raw->insData.gt,raw->insData.gtPre))
    {
        trace(3,"IMU data: data duplication --- %.4f.\n", raw->insData.gt);
        return 0;
    }
    else
    {    /* current and precious time difference is too large */
        trace(3, "IMU data: previous %.4f is greater than current %.4f.\n", raw->insData.gtPre, raw->insData.gt);
        return 0;
    }
    /* add imu measurement data */
    addimudata(&raw->imut,&raw->imu);
    return raw->insData.gt>0.0?4:0;

}


/* input imu raw data in backward--------------------------------------------*/
static int nextimub(const imu_t *imu,imu_t *data,int *index)
{
    int i; data->n=0;
    for (i=0;i<100;i++) addimudata(data,&imu->data[(*index)--]);
    if (data->n) return 4; else return 0;
}
/* input imu raw data for backward solution----------------------------------*/
static int decode_imu_m39b(raw_t *raw, unsigned char data)
{
    if (raw->imub.n==0) {
        prcopt_t *opt=(prcopt_t *)raw->optp;
        stream_t *str=(stream_t *)raw->strp;
        readimub(str->path,&raw->imub,opt->insopt.imudecfmt,opt->insopt.imuformat,
                 opt->insopt.imucoors,opt->insopt.imuvalfmt);

        raw->curb=raw->imub.n-1;
    }
    raw->imut.n=0;
    return nextimub(&raw->imub,&raw->imut,&raw->curb);
}
static int decode_imu_zhufengb(raw_t *raw, unsigned char data)
{
    if (raw->imub.n==0) {
        prcopt_t *opt=(prcopt_t *)raw->optp;
        stream_t *str=(stream_t *)raw->strp;
        readimub(str->path,&raw->imub,opt->insopt.imudecfmt,opt->insopt.imuformat,
                 opt->insopt.imucoors,opt->insopt.imuvalfmt);

        raw->curb=raw->imub.n-1;
    }
    raw->imut.n=0;
    return nextimub(&raw->imub,&raw->imut,&raw->curb);
}



/* input imu raw message -----------------------------------------------------
 * args   : raw_t *raw         IO     receiver raw data control struct
 *          unsigned char data I stream data (1 byte)
 * return >1:ok ,0:fail
 * --------------------------------------------------------------------------*/
extern int input_m39(raw_t *raw, unsigned char data)
{
    trace(3,"input_m39: data=%02X\n",data);
    raw->len=NUMBYTES_GI310;
    if (raw->dire) return decode_imu_m39b(raw,data);
    else           return decode_imu_m39 (raw,data);
}

extern int  input_zhufeng_raw(raw_t *raw, unsigned char data)
{
    trace(3," input_zhufeng_raw: data=%02X\n",data);
    raw->len=NUMBYTES_ZHUFENG;
    if (raw->dire) return decode_imu_zhufengb(raw,data);
    else           return decode_imu_zhufeng(raw,data);
}







/* input imu raw message from file --------------------------------------------
* input next imu raw message from file
* args   : raw_t  *raw   IO     receiver raw data control struct
*          FILE   *fp    I      file pointer
* return : status(-2: end of file, 1: ok 0: fail)
*----------------------------------------------------------------------------*/
extern int input_m39f(raw_t *raw,FILE *fp)
{
    int i,data,ret;

    trace(3,"input_imuf:\n");

    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_m39(raw,(unsigned char)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}
/* read imu measurement data from file-----------------------------------------
 * args   :  char *file     I  input imu measurement data file
 *           imu_t *imu     O  output imu measurement data
 *           int    decfmt  I  imu measurement data decode method
 *           int    imufmt  I  imu measurement data format
 *           int    coor    I  imu body coordinate frame
 *           int    valfmt  I  imu gyro measurement data format
 * return : number of imu measurement data
 * --------------------------------------------------------------------------*/
extern int readimub(const char *file,imu_t* imu,int decfmt,int imufmt,int coor,
                    int valfmt)
{
    FILE *fp;
    raw_t raw={0};
    int data,siz;

    trace(3,"readimub:\n");

    raw.imufmt=imufmt;
    imu->n=imu->nmax=0; imu->data=NULL;
    imu->format=decfmt; imu->coor=coor; imu->valfmt=valfmt;

    if (!(fp=fopen(file,"r"))) {
        trace(2,"imu measurement data file open error \n");
        return 0;
    }
    /* read imu message from file */
    while (1) {
        if ((data=fgetc(fp))==EOF) break;
        if ((input_m39(&raw,(unsigned char)data))) {

            if (imu->n>=imu->nmax) {
                trace(5,"readimub: imu->n=%d nmax=%d\n",imu->n,imu->nmax);

                siz=sizeof(imud_t)*(imu->nmax+=4096);
                if (!(imu->data=(imud_t *)realloc(imu->data,siz))) {

                    fprintf(stderr,"readimub :memory allocation error\n");
                    free(imu->data); imu->n=imu->nmax=0;
                    break;
                }
            }
            imu->data[imu->n++]=raw.imu;
        }
    }
    fclose(fp);
    return imu->n;
}
/* add imu measurement data -------------------------------------------------*/
extern int addimudata(imu_t *imu, const imud_t *data)
{
    imud_t *obs_data;

    if (imu->nmax<=imu->n) {
        if (imu->nmax<=0) imu->nmax=64; else imu->nmax*=2;
        if (!(obs_data=(imud_t *)realloc(imu->data,sizeof(imud_t)*imu->nmax))) {
            trace(1,"addimudata: memalloc error n=%dx%d\n",sizeof(imud_t),imu->nmax);
            free(imu->data); imu->data=NULL; imu->n=imu->nmax=0;
            return -1;
        }
        imu->data=obs_data;
    }
    imu->data[imu->n++]=*data;
    return 1;
}
/* free imu measurement data--------------------------------------------------*/
extern void freeimudata(imu_t *imu)
{
    trace(3,"freeimudata:\n");
    if (imu->data) {
        free(imu->data); imu->data=NULL; imu->n=imu->nmax=0;
    }
}
