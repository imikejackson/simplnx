#pragma once

#include <map>
#include <string>

// clang-format off
#include "Core/Filters/FindNeighborhoodsFilter.hpp"
#include "Core/Filters/RemoveFlaggedFeaturesFilter.hpp"
// @@__HEADER__TOKEN__DO__NOT__DELETE__@@

namespace complex
{
  static const std::map<complex::Uuid, complex::Uuid> k_SIMPL_to_Core
  {
    // syntax std::make_pair {Dream3d UUID , Dream3dnx UUID}, // dream3d-class-name
    {complex::Uuid::FromString("697ed3de-db33-5dd1-a64b-04fb71e7d63e").value(), complex::FilterTraits<FindNeighborhoodsFilter>::uuid}, // FindNeighborhoods
    {complex::Uuid::FromString("a8463056-3fa7-530b-847f-7f4cb78b8602").value(), complex::FilterTraits<RemoveFlaggedFeaturesFilter>::uuid}, // RemoveFlaggedFeatures
    // @@__MAP__UPDATE__TOKEN__DO__NOT__DELETE__@@
  };

} // namespace complex
// clang-format on
