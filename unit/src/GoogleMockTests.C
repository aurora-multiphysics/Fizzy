#include "gmock/gmock.h"

namespace {

class ValueSource {
public:
  virtual ~ValueSource() = default;
  virtual int value(int index) const = 0;
};

class MockValueSource : public ValueSource {
public:
  MOCK_METHOD(int, value, (int index), (const, override));
};

// Exercise GoogleMock's compiled implementation alongside MOOSE's GoogleTest.
TEST(GoogleMockIntegration, MatchesArgumentsAndReturnsConfiguredValue) {
  ::testing::StrictMock<MockValueSource> source;
  EXPECT_CALL(source, value(3)).Times(1).WillOnce(::testing::Return(42));

  const ValueSource &interface = source;
  EXPECT_EQ(interface.value(3), 42);
}

} // namespace
