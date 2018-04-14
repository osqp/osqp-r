OSQP_OUT_INCLUDE_DIR=./include/osqp/
OSQP_OUT_LIB_DIR=./lib/
OSQP_STATIC_LIB=$(OSQP_OUT_LIB_DIR)/libosqpstatic.a

OSQP_RCONFIG_DEFAULT_INCLUDE=../osqp_configure_R.h
OSQP_SRC_DIR=../osqp_sources/
OSQP_SRC_INCLUDE_DIR=$(OSQP_SRC_DIR)/include
OSQP_SRC_CONFIGURE_H=$(OSQP_SRC_INCLUDE_DIR)/osqp_configure.h

PARDISO_DIR=$(OSQP_SRC_DIR)/lin_sys/direct/pardiso/
SUITESPARSE_DIR=$(OSQP_SRC_DIR)/lin_sys/direct/suitesparse/

OSQP_SRC_INCLUDES_ALL= \
-I$(OSQP_SRC_INCLUDE_DIR) \
-I$(OSQP_SRC_DIR)/lin_sys \
-I$(PARDISO_DIR)      \
-I$(SUITESPARSE_DIR)  \
-I$(SUITESPARSE_DIR)/amd/include \
-I$(SUITESPARSE_DIR)/ldl/include

# We don't use wildcards since R complains about GNU make extensions.
# Otherwise we would do this:
# OSQP_SOURCES = $(wildcard $(OSQP_SRC_DIR)/src/*.c $(OSQP_SRC_DIR)/lin_sys/*.c $(PARDISO_DIR)/*.c $(SUITESPARSE_DIR)/*.c $(SUITESPARSE_DIR)/amd/src/*.c $(SUITESPARSE_DIR)/ldl/src/*.c)

OSQP_SOURCES = \
$(OSQP_SRC_DIR)/src/auxil.c   \
$(OSQP_SRC_DIR)/src/cs.c      \
$(OSQP_SRC_DIR)/src/ctrlc.c   \
$(OSQP_SRC_DIR)/src/kkt.c     \
$(OSQP_SRC_DIR)/src/lin_alg.c \
$(OSQP_SRC_DIR)/src/lin_sys.c \
$(OSQP_SRC_DIR)/src/osqp.c    \
$(OSQP_SRC_DIR)/src/polish.c  \
$(OSQP_SRC_DIR)/src/proj.c    \
$(OSQP_SRC_DIR)/src/scaling.c \
$(OSQP_SRC_DIR)/src/util.c    \
$(OSQP_SRC_DIR)/lin_sys/lib_handler.c         \
$(PARDISO_DIR)/pardiso.c   \
$(PARDISO_DIR)/pardiso_loader.c           \
$(SUITESPARSE_DIR)/SuiteSparse_config.c   \
$(SUITESPARSE_DIR)/suitesparse_ldl.c      \
$(SUITESPARSE_DIR)/ldl/src/ldl.c          \
$(SUITESPARSE_DIR)/amd/src/amd_1.c            \
$(SUITESPARSE_DIR)/amd/src/amd_2.c            \
$(SUITESPARSE_DIR)/amd/src/amd_aat.c          \
$(SUITESPARSE_DIR)/amd/src/amd_control.c      \
$(SUITESPARSE_DIR)/amd/src/amd_defaults.c     \
$(SUITESPARSE_DIR)/amd/src/amd_info.c         \
$(SUITESPARSE_DIR)/amd/src/amd_order.c        \
$(SUITESPARSE_DIR)/amd/src/amd_post_tree.c    \
$(SUITESPARSE_DIR)/amd/src/amd_postorder.c    \
$(SUITESPARSE_DIR)/amd/src/amd_preprocess.c   \
$(SUITESPARSE_DIR)/amd/src/amd_valid.c

OSQP_OBJECTS = $(OSQP_SOURCES:.c=.o)

CXX_STD=CXX11

all: $(OSQP_STATIC_LIB) osqp_includes

$(OSQP_OBJECTS): $(OSQP_SRC_CONFIGURE_H)
	$(CC) $(OSQP_SRC_INCLUDES_ALL) -I$(R_INCLUDE_DIR) $(CPICFLAGS) $(OSQP_FLAGS) $(CFLAGS) -c $(@:.o=.c) -o $@


$(OSQP_STATIC_LIB): $(OSQP_OBJECTS)
	mkdir -p $(OSQP_OUT_LIB_DIR)
	$(AR) $(ARFLAGS) $(OSQP_STATIC_LIB) $(OSQP_OBJECTS)
	$(RANLIB) $(OSQP_STATIC_LIB)

osqp_includes:
	mkdir -p $(OSQP_OUT_INCLUDE_DIR)
	cp $(OSQP_SRC_INCLUDE_DIR)/*.h $(OSQP_OUT_INCLUDE_DIR)

$(OSQP_SRC_CONFIGURE_H):
	cp $(OSQP_RCONFIG_DEFAULT_INCLUDE) $(OSQP_SRC_CONFIGURE_H)

clean:
	$(RM) $(OSQP_OBJECTS)
