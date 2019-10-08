#ifndef NEML_INTERFACE_H
#define NEML_INTERFACE_H

#include "nemlerror.h"
#include "models.h"
#include "nemlmath.h"
#include "singlecrystal.h"
#include "crystallography.h"
#include "rotations.h"


namespace neml {
std::unique_ptr<NEMLModel> parse_xml_unique(std::string fname, std::string mname);
}

#endif // NEML_INTERFACE_H
