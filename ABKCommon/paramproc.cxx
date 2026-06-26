/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, igor.markov1@gmail.com
***  Original Affiliation:   UCLA, Computer Science Department,
***                          Los Angeles, CA 90095-1596 USA
***
***  Permission is hereby granted, free of charge, to any person obtaining
***  a copy of this software and associated documentation files (the
***  "Software"), to deal in the Software without restriction, including
***  without limitation
***  the rights to use, copy, modify, merge, publish, distribute, sublicense,
***  and/or sell copies of the Software, and to permit persons to whom the
***  Software is furnished to do so, subject to the following conditions:
***
***  The above copyright notice and this permission notice shall be included
***  in all copies or substantial portions of the Software.
***
*** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***
***************************************************************************/

// Created: June 15, 1997 Igor Markov

// 970622 ilm    added NOPARAM
// 970820 ilm    moved enum ParamType inside class Param
//               added new ctor for use with NOPARAM
//               allowed for +option as well as -option
//               added member bool on() to tell if it was +option
//               added abkfatal(found(),...) when reading param, value
// 970824 ilm    changed the uninitialized value for int to MAX_INT
// 980313 ilm    fixed const-correctness

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "paramproc.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <strings.h>

#include <string>

#include "abkassert.h"
#include "abkstring.h"
static char _uninitialized[] = "Uninitialized";

namespace {

void failParamParse(const char* key, const char* value, const char* typeName) {
  std::string msg = " Invalid ";
  msg += typeName;
  msg += " value for command line parameter -";
  msg += key;
  msg += ": ";
  msg += value ? value : "<missing>";
  msg += "\n";
  abkfatal(0, msg.c_str());
}

void requireParamValue(const char* key, const char* value, const char* typeName) {
  if (value == NULL || value[0] == '\0') failParamParse(key, value, typeName);
}

bool isDecimalDigit(char c) {
  return isdigit(static_cast<unsigned char>(c)) != 0;
}

bool tokenLooksLikeAnotherOption(Param::Type type, const char* value) {
  if (value == NULL || value[0] == '\0') return false;
  if (value[0] != '-' && value[0] != '+') return false;

  const char second = value[1];
  switch (type) {
    case Param::INT:
      return !isDecimalDigit(second);
    case Param::UNSIGNED:
      return !isDecimalDigit(second);
    case Param::DOUBLE:
      return !isDecimalDigit(second) && second != '.';
    default:
      return false;
  }
}

int parseIntParam(const char* key, const char* value) {
  requireParamValue(key, value, "integer");
  errno = 0;
  char* end = NULL;
  const long parsed = strtol(value, &end, 10);
  if (end == value || *end != '\0' || errno == ERANGE || parsed < INT_MIN ||
      parsed > INT_MAX)
    failParamParse(key, value, "integer");
  return static_cast<int>(parsed);
}

unsigned parseUnsignedParam(const char* key, const char* value) {
  requireParamValue(key, value, "unsigned integer");
  errno = 0;
  char* end = NULL;
  const unsigned long parsed = strtoul(value, &end, 10);
  if (value[0] == '-' || end == value || *end != '\0' || errno == ERANGE ||
      parsed > UINT_MAX)
    failParamParse(key, value, "unsigned integer");
  return static_cast<unsigned>(parsed);
}

double parseDoubleParam(const char* key, const char* value) {
  requireParamValue(key, value, "floating point");
  errno = 0;
  char* end = NULL;
  const double parsed = strtod(value, &end);
  if (end == value || *end != '\0' || errno == ERANGE)
    failParamParse(key, value, "floating point");
  return parsed;
}

}  // namespace

Param::Param(Type pt, int argc, const char* const argv[])
    : _b(false),
      _on(false),
      _i(INT_MAX),
      _u(unsigned(-1)),
      _d(-1.29384756657),
      _s(_uninitialized),
      _pt(pt),
      _key("") {
  abkfatal(pt == NOPARAM,
           " This constructor can only work with Param:NOPARAM\n");
  (void)argv;  // please compiler

  _b = (argc < 2 ? true : false);
  return;
}

Param::Param(const char* key, Type pt, int argc, const char* const argv[])
    : _b(false),
      _on(false),
      _i(-1),
      _u(unsigned(-1)),
      _d(-1.23456),
      _s(_uninitialized),
      _pt(pt),
      _key(key) {
  abkfatal(strlen(_key) > 0, " Zero length key for command line parameter");

  int n = 0;
  if (_pt == NOPARAM) {
    if (argc < 2)
      _b = true;
    else
      _b = false;
    return;
  }
  while (++n < argc && !found()) {
    if (argv[n][0] == '-' || argv[n][0] == '+') {
      const char* start = argv[n] + 1;
      if (argv[n][0] == '-') {
        if (argv[n][1] == '-') start++;
      } else
        _on = true;

      if (strcasecmp(start, _key) == 0) {
        _b = true;
        if (n + 1 < argc) {
          switch (_pt) {
            case BOOL:
              break;
            case INT:
              if (tokenLooksLikeAnotherOption(_pt, argv[n + 1]))
                _i = 0;
              else
                _i = parseIntParam(_key, argv[n + 1]);
              break;
            case UNSIGNED:
              if (tokenLooksLikeAnotherOption(_pt, argv[n + 1]))
                _u = 0;
              else
                _u = parseUnsignedParam(_key, argv[n + 1]);
              break;
            case DOUBLE:
              if (tokenLooksLikeAnotherOption(_pt, argv[n + 1]))
                _d = 0.0;
              else
                _d = parseDoubleParam(_key, argv[n + 1]);
              break;
            case STRING:
              _s = argv[n + 1];
              break;
            default:
              abkfatal(0, " Unknown command line parameter");
          }
        }
      }
    }
  }
}

bool Param::found() const { return _b; }

bool Param::on() const  // true for +option, false otherwise
{
  abkfatal(found(), " Parameter not found: you need to check for this first\n");
  return _on;
}

int Param::getInt() const {
  abkfatal(_pt == INT, " Parameter is not INT ");
  abkfatal(found(), " Parameter not found: you need to check for this first\n");
  return _i;
}

unsigned Param::getUnsigned() const {
  abkfatal(_pt == UNSIGNED, " Parameter is not UNSIGNED ");
  abkfatal(found(),
           " UNSIGNED Parameter not found: you need to check for this first\n");
  return _u;
}

double Param::getDouble() const {
  abkfatal(_pt == DOUBLE, " Parameter is not DOUBLE ");
  abkfatal(found(),
           " DOUBLE parameter not found: you need to check for this first\n");
  return _d;
}

const char* Param::getString() const {
  abkfatal(_pt == STRING, " Parameter is not STRING");
  abkfatal(found(),
           " STRING parameter not found: you need to check for this first\n");
  return _s;
}
