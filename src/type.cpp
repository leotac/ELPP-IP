#include "type.h"

string which(ElppForm form)
{
   switch(form)
   {
      case NONE:
         return "None";
      case MCF:
         return "MCF";
      case SC:
         return "SC";
      case MinCut:
         return "MinCut";
      case SF:
         return "SF";
      case RLT:
         return "RLT";
      case MTZ:
         return "MTZ";
      case DL:
         return "DL";
      case DFJ:
         return "DFJ";
      case MCFsep:
         return "MCFsep";
      default:
         return "Unkn";
   }
}

