.SUFFIXES: .3 .3.html

CFLAGS		?= -g -W -Wall -Wstrict-prototypes -Wno-unused-parameter -Wwrite-strings
#If you're on Mac OSX without XQuarts, you'll need /usr/X11 instead of /opt/X11!
CPPFLAGS	+= `pkg-config --cflags --silence-errors cairo || echo '-I/opt/X11/include/cairo'`
VERSION		= 0.1.14
LDADD		= `pkg-config --libs --silence-errors cairo || echo '-L/opt/X11/lib -lcairo'`
#If you're on GNU/Linux, you'll need to uncomment this.
LDADD		+= `pkg-config --libs libbsd`
VERSIONS	= version_0_1_4.xml \
		  version_0_1_5.xml \
		  version_0_1_6.xml \
		  version_0_1_7.xml \
		  version_0_1_8.xml \
		  version_0_1_9.xml \
		  version_0_1_10.xml \
		  version_0_1_11.xml \
		  version_0_1_12.xml \
		  version_0_1_13.xml \
		  version_0_1_14.xml
EXAMPLES	= example0 \
		  example1 \
		  example2 \
		  example3 \
		  example4 \
		  example5 \
		  example6 \
		  example7 \
		  example8 \
		  example9 \
		  example10 \
		  example11 \
		  example12
PNGS		= example0.png \
		  example1.png \
		  example2.png \
		  example3.png \
		  example4.png \
		  example5.png \
		  example6.png \
		  example7.png \
		  example8.png \
		  example9.png \
		  example10.png \
		  example11.png \
		  example12.png
SRCS		= Makefile \
		  compat.post.h \
		  compat.pre.h \
		  compat.sh \
		  colours.c \
		  extern.h \
		  kplot.h \
		  array.c \
		  border.c \
		  bucket.c \
		  buffer.c \
		  draw.c \
		  example0.c \
		  example1.c \
		  example2.c \
		  example3.c \
		  example4.c \
		  example5.c \
		  example6.c \
		  example7.c \
		  example8.c \
		  example9.c \
		  example10.c \
		  example11.c \
		  example12.c \
		  grid.c \
		  hist.c \
		  label.c \
		  kdata.c \
		  kplot.c \
		  margin.c \
		  mean.c \
		  plotctx.c \
		  reallocarray.c \
		  stddev.c \
		  test-reallocarray.c \
		  tic.c \
		  vector.c
OBJS		= array.o \
		  border.o \
		  bucket.o \
		  buffer.o \
		  colours.o \
		  draw.o \
		  grid.o \
		  hist.o \
		  label.o \
		  kdata.o \
		  kplot.o \
		  margin.o \
		  mean.o \
		  plotctx.o \
		  reallocarray.o \
		  stddev.o \
		  tic.o \
		  vector.o
PREFIX		= /usr
HTMLS		= index.html \
		  man/kdata_array_alloc.3.html \
		  man/kdata_array_fill.3.html \
		  man/kdata_bucket_alloc.3.html \
		  man/kdata_bucket_add.3.html \
		  man/kdata_buffer_alloc.3.html \
		  man/kdata_buffer_copy.3.html \
		  man/kdata_destroy.3.html \
		  man/kdata_get.3.html \
		  man/kdata_hist_alloc.3.html \
		  man/kdata_hist_add.3.html \
		  man/kdata_mean_alloc.3.html \
		  man/kdata_mean_attach.3.html \
		  man/kdata_pmfmean.3.html \
		  man/kdata_pmfstddev.3.html \
		  man/kdata_pmfvar.3.html \
		  man/kdata_stddev_alloc.3.html \
		  man/kdata_stddev_attach.3.html \
		  man/kdata_vector_append.3.html \
		  man/kdata_vector_alloc.3.html \
		  man/kdata_xmax.3.html \
		  man/kdata_xmean.3.html \
		  man/kdata_xmin.3.html \
		  man/kdata_xstddev.3.html \
		  man/kdatacfg_defaults.3.html \
		  man/kplot.3.html \
		  man/kplot_alloc.3.html \
		  man/kplot_attach_data.3.html \
		  man/kplot_attach_datas.3.html \
		  man/kplot_attach_smooth.3.html \
		  man/kplot_detach.3.html \
		  man/kplot_draw.3.html \
		  man/kplot_free.3.html \
		  man/kplot_get_datacfg.3.html \
		  man/kplot_get_plotcfg.3.html \
	 	  man/kplotcfg_defaults.3.html
MANS		= man/kdata_array_alloc.3 \
		  man/kdata_array_fill.3 \
		  man/kdata_bucket_alloc.3 \
		  man/kdata_bucket_add.3 \
		  man/kdata_buffer_alloc.3 \
		  man/kdata_buffer_copy.3 \
		  man/kdata_destroy.3 \
		  man/kdata_get.3 \
		  man/kdata_hist_alloc.3 \
		  man/kdata_hist_add.3 \
		  man/kdata_mean_alloc.3 \
		  man/kdata_mean_attach.3 \
		  man/kdata_pmfmean.3 \
		  man/kdata_pmfstddev.3 \
		  man/kdata_pmfvar.3 \
		  man/kdata_stddev_alloc.3 \
		  man/kdata_stddev_attach.3 \
		  man/kdata_vector_append.3 \
		  man/kdata_vector_alloc.3 \
		  man/kdata_xmax.3 \
		  man/kdata_xmean.3 \
		  man/kdata_xmin.3 \
		  man/kdata_xstddev.3 \
		  man/kdatacfg_defaults.3 \
		  man/kplot.3 \
		  man/kplot_alloc.3 \
		  man/kplot_attach_data.3 \
		  man/kplot_attach_datas.3 \
		  man/kplot_attach_smooth.3 \
		  man/kplot_detach.3 \
		  man/kplot_draw.3 \
		  man/kplot_free.3 \
		  man/kplot_get_datacfg.3 \
		  man/kplot_get_plotcfg.3 \
	 	  man/kplotcfg_defaults.3

all: libkplot.a

install: all

uninstall: ;

www: $(HTMLS) $(PNGS) kplot.tgz kplot.tgz.sha512

installwww: www
	mkdir -p $(PREFIX)
	mkdir -p $(PREFIX)/snapshots
	install -m 0444 $(PNGS) $(HTMLS) index.css $(PREFIX)
	install -m 0444 kplot.tgz $(PREFIX)/snapshots/kplot-$(VERSION).tgz
	install -m 0444 kplot.tgz.sha512 $(PREFIX)/snapshots/kplot-$(VERSION).tgz.sha512
	install -m 0444 kplot.tgz $(PREFIX)/snapshots
	install -m 0444 kplot.tgz.sha512 $(PREFIX)/snapshots

$(EXAMPLES): libkplot.a

EXAMPLE_LIBS = libkplot.a $(LDADD) -lm 

example0: example0.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example1: example1.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example2: example2.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example3: example3.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example4: example4.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example5: example5.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example6: example6.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example7: example7.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example8: example8.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example9: example9.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example10: example10.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example11: example11.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example12: example12.c libkplot.a
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(EXAMPLE_LIBS)

example0.png: example0
	./example0

example1.png: example1
	./example1

example2.png: example2
	./example2

example3.png: example3
	./example3

example4.png: example4
	./example4

example5.png: example5
	./example5

example6.png: example6
	./example6

example7.png: example7
	./example7

example8.png: example8
	./example8

example9.png: example9
	./example9

example10.png: example10
	./example10

example11.png: example11
	./example11

example12.png: example12
	./example12

libkplot.a: $(OBJS)
	$(AR) rs $@ $(OBJS)

$(OBJS): kplot.h compat.h extern.h

compat.h: compat.pre.h compat.post.h
	( cat compat.pre.h ; \
	  CC="$(CC)" CFLAGS="$(CFLAGS)" sh compat.sh reallocarray ; \
	  cat compat.post.h ; ) >$@

.3.3.html:
	mandoc -Thtml -Oman=%N.%S.html $< >$@

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

index.html: index.xml $(VERSIONS)
	sblg -t index.xml -o- $(VERSIONS) | sed "s!@VERSION@!$(VERSION)!g" >$@

kplot.tgz.sha512: kplot.tgz
	openssl dgst -sha512 kplot.tgz >$@

kplot.tgz:
	mkdir -p .dist/kplot-$(VERSION)
	mkdir -p .dist/kplot-$(VERSION)/man
	cp $(SRCS) .dist/kplot-$(VERSION)
	cp $(MANS) .dist/kplot-$(VERSION)/man
	(cd .dist && tar zcf ../$@ kplot-$(VERSION))
	rm -rf .dist

clean:
	rm -f libkplot.a compat.h test-reallocarray 
	rm -f $(EXAMPLES)
	rm -rf *.dSYM
	rm -f $(OBJS)
	rm -f $(HTMLS) $(PNGS)
	rm -f kplot.tgz kplot.tgz.sha512
