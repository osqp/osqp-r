OSQP_OUT_INCLUDE_DIR=../include/osqp/
OSQP_OUT_LIB_DIR=../lib/
OSQP_STATIC_LIB=$(OSQP_OUT_LIB_DIR)/libosqp.a

OSQP_RCONFIG_DEFAULT_INCLUDE=../osqp_configure_R.h
QDLDL_RCONFIG_DEFAULT_INCLUDE=../qdldl_types_R.h
PRINTING_R_INCLUDE=../printing_R.h

OSQP_SRC_DIR=../osqp_sources
OSQP_SRC_CONFIGURE_H=$(OSQP_SRC_DIR)/include/private/osqp_configure.h
OSQP_SRC_PRINTING_H=$(OSQP_SRC_DIR)/include/private/printing.h

QDLDL_VENDORED_DIR=../qdldl_sources
QDLDL_SRC_TYPES_H=$(QDLDL_VENDORED_DIR)/include/qdldl_types.h

ALGEBRA_COMMON_DIR=$(OSQP_SRC_DIR)/algebra/_common
QDLDL_IFACE_DIR=$(ALGEBRA_COMMON_DIR)/lin_sys/qdldl
AMD_SRC_DIR=$(QDLDL_IFACE_DIR)/amd/src
ALGEBRA_BUILTIN_DIR=$(OSQP_SRC_DIR)/algebra/builtin

OSQP_SRC_INCLUDES_ALL= \
-I$(OSQP_SRC_DIR)/include/public \
-I$(OSQP_SRC_DIR)/include/private \
-I$(ALGEBRA_COMMON_DIR) \
-I$(QDLDL_IFACE_DIR) \
-I$(QDLDL_IFACE_DIR)/amd/include \
-I$(ALGEBRA_BUILTIN_DIR) \
-I$(QDLDL_VENDORED_DIR)/include

# OSQP core sources
OSQP_SOURCES = \
$(OSQP_SRC_DIR)/src/osqp_api.c \
$(OSQP_SRC_DIR)/src/auxil.c \
$(OSQP_SRC_DIR)/src/error.c \
$(OSQP_SRC_DIR)/src/polish.c \
$(OSQP_SRC_DIR)/src/scaling.c \
$(OSQP_SRC_DIR)/src/util.c \
$(OSQP_SRC_DIR)/src/profilers.c

# Platform-specific timing
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  OSQP_SOURCES += $(OSQP_SRC_DIR)/src/timing_macos.c
  OSQP_SOURCES += $(OSQP_SRC_DIR)/src/interrupt_unix.c
else ifeq ($(UNAME_S),Linux)
  OSQP_SOURCES += $(OSQP_SRC_DIR)/src/timing_linux.c
  OSQP_SOURCES += $(OSQP_SRC_DIR)/src/interrupt_unix.c
else
  OSQP_SOURCES += $(OSQP_SRC_DIR)/src/timing_windows.c
  OSQP_SOURCES += $(OSQP_SRC_DIR)/src/interrupt_windows.c
endif

# Algebra common sources
OSQP_SOURCES += \
$(ALGEBRA_COMMON_DIR)/csc_math.c \
$(ALGEBRA_COMMON_DIR)/csc_utils.c \
$(ALGEBRA_COMMON_DIR)/kkt.c \
$(ALGEBRA_COMMON_DIR)/reduced_kkt.c

# QDLDL interface
OSQP_SOURCES += $(QDLDL_IFACE_DIR)/qdldl_interface.c

# AMD sources
OSQP_SOURCES += \
$(AMD_SRC_DIR)/SuiteSparse_config.c \
$(AMD_SRC_DIR)/amd_1.c \
$(AMD_SRC_DIR)/amd_2.c \
$(AMD_SRC_DIR)/amd_aat.c \
$(AMD_SRC_DIR)/amd_control.c \
$(AMD_SRC_DIR)/amd_defaults.c \
$(AMD_SRC_DIR)/amd_info.c \
$(AMD_SRC_DIR)/amd_order.c \
$(AMD_SRC_DIR)/amd_post_tree.c \
$(AMD_SRC_DIR)/amd_postorder.c \
$(AMD_SRC_DIR)/amd_preprocess.c \
$(AMD_SRC_DIR)/amd_valid.c

# Builtin algebra sources
OSQP_SOURCES += \
$(ALGEBRA_BUILTIN_DIR)/algebra_libs.c \
$(ALGEBRA_BUILTIN_DIR)/matrix.c \
$(ALGEBRA_BUILTIN_DIR)/vector.c

# Vendored QDLDL source
OSQP_SOURCES += $(QDLDL_VENDORED_DIR)/src/qdldl.c

OSQP_OBJECTS = $(OSQP_SOURCES:.c=.o)

all: $(OSQP_STATIC_LIB) osqp_includes

$(OSQP_OBJECTS): copy_r_headers
	$(CC) -c $(@:.o=.c) -o $@ $(OSQP_SRC_INCLUDES_ALL) -I$(R_INCLUDE_DIR) $(CPICFLAGS) $(OSQP_FLAGS) $(CFLAGS)

.PHONY: copy_r_headers
copy_r_headers:
	cp $(OSQP_RCONFIG_DEFAULT_INCLUDE) $(OSQP_SRC_CONFIGURE_H)
	cp $(QDLDL_RCONFIG_DEFAULT_INCLUDE) $(QDLDL_SRC_TYPES_H)
	cp $(PRINTING_R_INCLUDE) $(OSQP_SRC_PRINTING_H)

$(OSQP_STATIC_LIB): $(OSQP_OBJECTS)
	mkdir -p $(OSQP_OUT_LIB_DIR)
	$(AR) $(ARFLAGS) $(OSQP_STATIC_LIB) $(OSQP_OBJECTS)
	$(RANLIB) $(OSQP_STATIC_LIB)

osqp_includes:
	mkdir -p $(OSQP_OUT_INCLUDE_DIR)
	cp $(OSQP_SRC_DIR)/include/public/*.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(OSQP_SRC_DIR)/include/private/*.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(ALGEBRA_BUILTIN_DIR)/algebra_impl.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(ALGEBRA_COMMON_DIR)/csc_math.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(ALGEBRA_COMMON_DIR)/csc_utils.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(ALGEBRA_COMMON_DIR)/kkt.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(ALGEBRA_COMMON_DIR)/reduced_kkt.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(QDLDL_IFACE_DIR)/qdldl_interface.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(QDLDL_VENDORED_DIR)/include/qdldl.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(QDLDL_VENDORED_DIR)/include/qdldl_types.h $(OSQP_OUT_INCLUDE_DIR)
	cp $(QDLDL_VENDORED_DIR)/include/qdldl_version.h $(OSQP_OUT_INCLUDE_DIR)


clean:
	$(RM) $(OSQP_OBJECTS)
