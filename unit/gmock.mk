# Build GoogleMock against MOOSE's GoogleTest, without compiling a second gtest.
GMOCK_DIR := $(abspath $(CURRENT_DIR)/../contrib/googletest/googlemock)
GMOCK_BUILD_DIR := $(CURRENT_DIR)/build/gmock/$(libmesh_HOST).$(METHOD)
GMOCK_OBJECT := $(GMOCK_BUILD_DIR)/gmock-all.lo
GMOCK_LIB := $(GMOCK_BUILD_DIR)/libgmock.la

# Cleaning must still work before submodules have been initialized.
ifeq ($(wildcard $(GMOCK_DIR)/src/gmock-all.cc),)
ifneq ($(filter-out clean clobber cleanall clobberall,$(or $(MAKECMDGOALS),all)),)
$(error GoogleMock is missing. Run 'git submodule update --init contrib/googletest' from the Fizzy repository root)
endif
endif

ADDITIONAL_INCLUDES += -I$(GMOCK_DIR)/include
ADDITIONAL_LIBS += $(GMOCK_LIB)
ADDITIONAL_DEPEND_LIBS += $(GMOCK_LIB)

$(GMOCK_BUILD_DIR):
	mkdir -p "$@"

$(GMOCK_OBJECT): $(GMOCK_DIR)/src/gmock-all.cc $(CURRENT_DIR)/gmock.mk | $(GMOCK_BUILD_DIR)
	$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=compile --quiet \
	  $(libmesh_CXX) $(libmesh_CPPFLAGS) $(CXXFLAGS) $(libmesh_CXXFLAGS) \
	  $(ADDITIONAL_CPPFLAGS) -I$(FRAMEWORK_DIR)/contrib/gtest \
	  -I$(GMOCK_DIR)/include -I$(GMOCK_DIR) \
	  -MMD -MP -MF $@.d -MT $@ -c $< -o $@

$(GMOCK_LIB): $(GMOCK_OBJECT) $(gtest_LIB)
	$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) -o $@ $(GMOCK_OBJECT) $(gtest_LIB) \
	  $(LDFLAGS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(GMOCK_BUILD_DIR)
	$(libmesh_LIBTOOL) --mode=install --quiet install -c $@ $(GMOCK_BUILD_DIR)

# MOOSE's clean/clobber targets remove unit/build, including these artifacts.
-include $(GMOCK_OBJECT).d
