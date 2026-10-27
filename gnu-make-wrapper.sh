#!/bin/sh

mode=$1
pkgdir=$2
shift 2

helper_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
host_system=${UMPACK_HOST_SYSTEM:-$(uname -s)}
case "$host_system" in
Darwin)
	cxx="$helper_dir/cxx-wrapper.sh"
	cc=clang
	;;
*)
	cxx=g++
	cc=gcc
	;;
esac

case "$mode" in
make-gnu)
	makefile="./make-gnu"
	optflags='-g -DABKDEBUG'
	ccflags='-fPIC -Wall -pedantic -pipe'
	;;
makeOpt-gnu)
	makefile="./makeOpt-gnu"
	optflags='-O3 -march=native -funroll-loops $(ABKSYMBOLS)'
	ccflags='-fPIC -pipe'
	;;
*)
	echo "gnu-make-wrapper.sh: unknown mode '$mode'" >&2
	exit 1
	;;
esac

make_pkg() {
	case "$1" in
	with-c)
		shift
		exec make \
			MAKE="$makefile" DEPENDOPT="-MM" CC4LIBS="$cxx" \
			OPTFLAGS="$optflags" CCFLAGS="$ccflags" \
			C="$cc \$(OPTFLAGS) \$(CCFLAGS)" \
			CC="$cxx \$(OPTFLAGS) \$(CCFLAGS)" \
			LD="$cxx -shared" LDFINAL="$cxx" SODIRKEY="-Wl,-rpath," \
			SOLABEL="" \
			"$@"
		;;
	*)
		exec make \
			MAKE="$makefile" DEPENDOPT="-MM" CC4LIBS="$cxx" \
			OPTFLAGS="$optflags" CCFLAGS="$ccflags" \
			CC="$cxx \$(OPTFLAGS) \$(CCFLAGS)" \
			LD="$cxx -shared" LDFINAL="$cxx" SODIRKEY="-Wl,-rpath," \
			SOLABEL="" \
			"$@"
		;;
	esac
}

make_pkg "$@"
