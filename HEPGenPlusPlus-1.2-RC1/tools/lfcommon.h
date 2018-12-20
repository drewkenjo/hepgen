#ifndef lfcommon_h
#define lfcommon_h

#define LFWR_SUCCESS 0
#define LFWR_ERROR -1
#define LFWR_EOF -2



namespace lfr {
/*! \struct lHeader
 *  \brief Holds the information of the header-part of the lepto-file
 */
typedef struct
{
    int lst[40];
    float parl[30];
    float cutl[14];
    float unidentified[20];
} lHeader;

/*! \struct lBeamParticle
 *  \brief holds the information of one lujet
 *  This is a lujet-type struct, containing k,p and lu2kine.
 */
typedef struct
{
    int k[5];
    float p[5];
    int lu2kine;
} lBeamParticle;

/*! \struct lEvent
 *  \brief This holds the whole Lepto-Event information
 *  in beamParts there will be an array of lBeamParticle depending on the size of the event.
 *  Make sure to free it after use!
 */
typedef struct
{
    int size;
    float uservar[20];
    int lst[40];
    float parl[30];
    float cut[14];
    float q2,w2,nu,x_bj,y;
    unsigned int nBeamParticle;
    lBeamParticle* beamParts;
} lEvent;

}

#endif
