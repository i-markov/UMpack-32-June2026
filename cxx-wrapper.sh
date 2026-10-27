#!/usr/bin/perl

use strict;
use warnings;

# The package Makefiles use GNU compiler and linker spellings. Keep those
# command lines unchanged on Linux, and translate the small incompatible
# subset when Apple Clang is the platform C++ driver.
my $host_system = $ENV{UMPACK_HOST_SYSTEM};
if ( !defined $host_system || $host_system eq "" ) {
    chomp( $host_system = `uname -s` );
}
exec "g++", @ARGV if $host_system ne "Darwin";

my @args;
foreach my $arg (@ARGV) {
    if ( $arg eq "-shared" ) {
        # GNU shared objects allow references that their eventual executable or
        # another library will provide. Match that behavior for legacy packages.
        push @args, "-dynamiclib", "-Wl,-undefined,dynamic_lookup";
    }
    elsif ( $arg eq "-static"
        || $arg eq "-Wl,--no-as-needed"
        || $arg eq "-Wl,--start-group"
        || $arg eq "-Wl,--end-group" )
    {
        next;
    }
    elsif ( $arg =~ /^-Wl,-soname,(.*)$/ ) {
        push @args, "-Wl,-install_name,\@rpath/$1";
    }
    else {
        push @args, $arg;
    }
}

exec "clang++", @args;
