%module lhepgen

%include std_string.i
using std::string;

%include "cpointer.i"


%pointer_class(int, intp);
%pointer_class(double, doublep);


%{
 /* Put header files here or function declarations like below */
 #include "hWeightingInterface.h"
 #include "reweightKine.h"
 #include "hphysicsgen.h"
 #include "hphigen.h"
 #include "hpigen.h"
 #include "homegagen.h"
 #include "hrhogen.h"
 #include "hrhoplusgen.h"
 #include "hMosseGen.h"
 #include "hdvcsgen.h"
 #include "myTHEO.hh"
 #include "libGKPi0.h"
 #include "gkSubProcessTable.h"
 #include "hvector.h"
 #include "hlorentzvector.h"
 #include "hrotmat.h"
 #include "hpamgen.h"
 #include "hjpsigen.h"
 %}


 %include "hvector.h"
 %include "hlorentzvector.h"
 %include "hWeightingInterface.h"
 %include "reweightKine.h"
 %include "hphysicsgen.h"
 %include "hphigen.h"
 %include "hpigen.h"
 %include "hrotmat.h"
 %include "hpionicdata.h"
 %include "homegagen.h"
 %include "hrhogen.h"
 %include "hrhoplusgen.h"
 %include "hjpsigen.h"
 %include "hpamgen.h"
 %include "hdvcsgen.h"
 %include "libGKPi0.h"
 %include "gkSubProcessTable.h"
%pointer_class(HPionicData, HPionicDataPointer);
%pointer_class(hWeightInterface, hWeightInterfacePointer);
 
