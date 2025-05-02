#ifndef TYPES_H
#define TYPES_H

#include "tensor.H"
#include "vector.H"
#include "quaternion.H"
#include "dictionary.H"
#include <string>

namespace sdfibm {

using Foam::label;
using Foam::scalar;
using Foam::vector;
using Foam::tensor;
using Foam::quaternion;
using Foam::dictionary;

constexpr scalar SMALL = 1e-6;

} // namespace sdfibm
#endif // TYPES_H
