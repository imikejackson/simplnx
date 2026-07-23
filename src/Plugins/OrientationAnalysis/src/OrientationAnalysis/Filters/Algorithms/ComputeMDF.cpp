#include "ComputeMDF.hpp"

using namespace nx::core;

// -----------------------------------------------------------------------------
ComputeMDF::ComputeMDF(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, ComputeMDFInputValues* inputValues)
: m_DataStructure(dataStructure)
, m_InputValues(inputValues)
, m_ShouldCancel(shouldCancel)
, m_MessageHandler(mesgHandler)
{
}

// -----------------------------------------------------------------------------
ComputeMDF::~ComputeMDF() noexcept = default;

// -----------------------------------------------------------------------------
Result<> ComputeMDF::operator()()
{
  // STUB: The KDE-based MDF computation is implemented in a subsequent task.
  return {};
}
