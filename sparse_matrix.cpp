
#include "sparse_matrix.h"


#ifdef RV_SUITESPARSE_SUPPORT
namespace RobotVision{
 CholmodSingleton CholmodSingleton::_instance;


  CholmodSingleton & CholmodSingleton::getInstance()
  {
    return _instance;
  }

}
#endif


