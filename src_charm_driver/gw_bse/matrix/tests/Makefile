include ../../config.mk

C_FILES=$(wildcard *.C)
CI_FILES=$(wildcard *.ci)
MODULES=$(CI_FILES:.ci=.decl.h)
OBJS=$(addprefix $(BUILDDIR)/, $(C_FILES:.C=.o))
EXES=$(BUILDDIR)/matrix_test

default: all
all: $(BUILDDIR) $(EXES)

$(BUILDDIR):
	mkdir $(BUILDDIR)

$(EXES): $(OBJS) $(BUILDDIR)/matrix.o
	$(CHARMC) $(OBJS) $(BUILDDIR)/matrix.o -o $@

make_ci: $(MODULES)

%.decl.h: %.ci
	$(CHARMC) $<

$(BUILDDIR)/matrix.o: ../matrix.C ../matrix.h ../matrix.ci
	make -C ../

$(BUILDDIR)/%.o: %.C %.h $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.C $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf charmrun
	rm -rf *.decl.h *.def.h
	rm -rf $(OBJS)
	rm -rf $(EXES)
