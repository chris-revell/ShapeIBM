//
//  LocalToGlobal.hpp
//  ImmersedBoundary
//
//  Created by Christopher Revell on 15/02/2019.
//
//

#ifndef LocalToGlobal_hpp
#define LocalToGlobal_hpp

#include <stdio.h>
#include "LocalToGlobal.hpp"
#include "element.hpp"
#include "cell.hpp"
#include "tissue.hpp"
#include <armadillo>
#include <vector>

void LocalToGlobal(tissue& Tissue);

#endif /* LocalToGlobal_hpp */
