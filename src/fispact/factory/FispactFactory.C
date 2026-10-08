#ifdef FIZZY_UNIT_TEST
#include "FispactContextMock.h"
#else
#include "DummyFispactContext.h"
#include "FispactContext.h"
#endif

#include "FispactFactory.h"

std::unique_ptr<FispactContextBase> createFispactContext(bool dummy_photons) {
#ifdef FIZZY_UNIT_TEST
  return std::make_unique<FispactContextMock>();
#else
  if (!dummy_photons)
    return std::make_unique<FispactContext>();
  else
    return std::make_unique<DummyFispactContext>();
#endif
}
