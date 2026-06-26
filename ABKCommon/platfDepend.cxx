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

// August   22, 1997   created   by Igor Markov  VLSI CAD UCLA ABKGROUP
// November 30, 1997   additions by Max Moroz,   VLSI CAD UCLA ABKGROUP

// note: not thread-safe

// This file to be included into all projects in the group
// it contains platform-specific code portability of which relies
// on symbols defined by compilers

// 970825 mro made corrections to conditional compilation in ctors for
//            Platform and User:
//            i)   _MSC_VER not __MSC_VER (starts with single underscore)
//            ii)  allocated more space for _infoLines
//            iii) changed nested #ifdefs to #if ... #elif ... #else
// 970923 ilm added abk_dump_stack()
// 971130 ilm added Max Moroz code for memory estimate and reworked
//            class MemUsage()

#include "abkassert.h"
#include "abktypes.h"
#include "infolines.h"

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#define _X86_
#include <windows.h>
#endif

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef linux
#include <sys/utsname.h>
#endif

#if defined(linux)
#include <pwd.h>
#include <stdlib.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <iomanip>
#include <iostream>
#endif

#ifdef _MSC_VER
#include <psapi.h>
#include <windows.h>
#endif

using std::ifstream;

#if defined(linux)
namespace {

bool isProcBlank(char c) { return c == ' ' || c == '\t'; }

bool isProcLineEnd(char c) { return c == '\0' || c == '\n' || c == '\r'; }

bool parseDoubleToken(const char* text, double* value) {
  if (text == NULL || value == NULL) return false;
  while (isProcBlank(*text)) ++text;
  if (*text == '\0') return false;

  errno = 0;
  char* end = NULL;
  const double parsed = strtod(text, &end);
  if (end == text || errno == ERANGE) return false;
  while (isProcBlank(*end)) ++end;
  if (*end != '\0') return false;

  *value = parsed;
  return true;
}

bool parseProcKilobytesLine(const char* line, const char* label, long* value) {
  if (line == NULL || label == NULL || value == NULL) return false;
  const size_t labelLen = strlen(label);
  if (strncmp(line, label, labelLen) != 0) return false;

  const char* number = line + labelLen;
  while (isProcBlank(*number)) ++number;
  errno = 0;
  char* end = NULL;
  const long parsed = strtol(number, &end, 10);
  if (end == number || errno == ERANGE || parsed < 0) return false;

  while (isProcBlank(*end)) ++end;
  if (end[0] == 'k' && end[1] == 'B') {
    end += 2;
    while (isProcBlank(*end)) ++end;
  }
  if (!isProcLineEnd(*end)) return false;

  *value = parsed;
  return true;
}

}  // namespace
#endif

/* ======================== IMPLEMENTATIONS ======================== */

#ifdef _MSC_VER

double PrintWindowsMemoryInfo(DWORD processID) {
  const unsigned kMegaByte = 1024 * 1024;
  HANDLE hProcess;
  PROCESS_MEMORY_COUNTERS pmc;

  // Print the process identifier.

  // printf( "\nProcess ID: %u\n", processID );

  // Print information about the memory usage of the process.

  hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE,
                         processID);
  if (NULL == hProcess) return -1.0;

  double memUsage = -1.0;
  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
    // printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
    // printf( "\tPeakWorkingSetSize: 0x%08X\n",
    //           pmc.PeakWorkingSetSize );
    // printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
    // printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n",
    //           pmc.QuotaPeakPagedPoolUsage );
    // printf( "\tQuotaPagedPoolUsage: 0x%08X\n",
    //           pmc.QuotaPagedPoolUsage );
    // printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n",
    //           pmc.QuotaPeakNonPagedPoolUsage );
    // printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n",
    //           pmc.QuotaNonPagedPoolUsage );
    // printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage );
    // printf( "\tPeakPagefileUsage: 0x%08X\n",
    //           pmc.PeakPagefileUsage );
    memUsage = pmc.WorkingSetSize;
  }

  CloseHandle(hProcess);
  return memUsage / kMegaByte;
}

double winMemUsage(void) {
  return PrintWindowsMemoryInfo(GetCurrentProcessId());
}

#endif

void abk_dump_stack() {
  return;
#if defined(linux)
  printf("\n  --- ABK GROUP DELIBERATE STACK DUMP  ----------- \n\n");
  fflush(stdout);
  char s[160];
  unsigned ProcId = getpid();
  snprintf(s, sizeof(s),
           "/bin/echo \"bt   \nquit   \\n\" | "
           "gdb -q /proc/%d/exe %d",
           ProcId, ProcId);
  int status = system(s);
  (void)status;
  printf("  ------------------ END STACK DUMP -------------- \n");
  fflush(stdout);
#else
// fprintf(stderr," abk_dump_stack(): Can't dump stack on this platform\n");
// fflush(stderr);
#endif
  return;
}

void abk_call_debugger() {
#if defined(linux)
  unsigned ProcId = getpid();
  printf("\n  --- ATTACHING DEBUGGER to process %d ", ProcId);
  printf(" (an ABKGROUP utility) --- \n\n");
  fflush(stdout);
  char s[160];
  snprintf(s, sizeof(s), "gdb -q /proc/%d/exe %d", ProcId, ProcId);
  int status = system(s);
  (void)status;
  printf("  ------------------ CONTINUING -------------- \n");
  fflush(stdout);
#else
  fprintf(stderr,
          " abk_call_debugger(): Can't call debugger on this platform\n");
  fflush(stderr);
#endif
  return;
}

double getPeakMemoryUsage() {
#if defined(_MSC_VER)
  return winMemUsage();
#elif defined(linux)
  char buf[1000];
  ifstream ifs("/proc/self/stat");
  for (unsigned k = 0; k != 23; k++) {
    if (!(ifs >> buf)) return -1.0;
  }
  //  cout << k << ": " << buf << endl; }
  double pages = -1.0;
  if (!parseDoubleToken(buf, &pages)) return -1.0;
  return (1.0 / (1024. * 1024.)) * pages;
#else
  return -1;
#endif
}

static double getMemoryUsageEstimate();

MemUsage::MemUsage() {
  _peak = getPeakMemoryUsage();
#if defined(linux)
  _estimate = _peak;
  return;
#endif
  _estimate = getMemoryUsageEstimate();
}

static double getVmRSS();
static double getTotalPhysMemSize();
static double getAvailablePhysMemSize();
static double getMinorPageFaults();
static double getMajorPageFaults();

double VMemUsage::_phys_total = getTotalPhysMemSize();
double VMemUsage::measureMinPageFaults() { return getMinorPageFaults(); }
double VMemUsage::measureMajPageFaults() { return getMajorPageFaults(); }
double VMemUsage::measureRSS() { return getVmRSS(); }
double VMemUsage::measurePhysAvailable() { return getAvailablePhysMemSize(); }

VMemUsage::VMemUsage() {
  _rss = getVmRSS();
  _minflt = getMinorPageFaults();
  _majflt = getMajorPageFaults();
  _phys = getAvailablePhysMemSize();
}

void MemChange::resetMark() { _memUsageMark = getMemoryUsageEstimate(); }

MemChange::MemChange() { resetMark(); }

MemChange::operator double() const {
#if defined(linux)
  return getPeakMemoryUsage() - _memUsageMark;
#else
  return -1.0;
#endif
}

Platform::Platform() {
#if defined(linux)
  struct utsname buf;
  uname(&buf);
  _infoLine = new char[strlen(buf.sysname) + strlen(buf.release) +
                       strlen(buf.version) + strlen(buf.machine) + 30];
  snprintf(_infoLine,
           strlen(buf.sysname) + strlen(buf.release) + strlen(buf.version) +
               strlen(buf.machine) + 30,
           "# Platform     : %s %s %s %s \n", buf.sysname, buf.release,
           buf.version, buf.machine);
#elif defined(_MSC_VER)
  _infoLine = new char[40];
  strcpy(_infoLine, "# Platform     : MS Windows \n");
#else
  _infoLine = new char[40];
  strcpy(_infoLine, "# Platform     : unknown \n");
#endif
}

#if !defined(linux) && !defined(_MSC_VER)
extern int gethostname(char *, int);
#endif

User::User() {
#if defined(linux)
  char host[100];
  gethostname(host, 31);
  struct passwd *pwr = getpwuid(getuid());
  if (pwr == NULL) {
    _infoLine = new char[40];
    strcpy(_infoLine, "# User         : unknown \n");
    return;
  }

  _infoLine = new char[strlen(pwr->pw_name) + 130];
  snprintf(_infoLine, strlen(pwr->pw_name) + 130, "# User         : %s@%s \n",
           pwr->pw_name, host);
#elif defined(_MSC_VER)
  _infoLine = new char[40];
  strcpy(_infoLine, "# User         : unknown \n");
#else
  _infoLine = new char[40];
  strcpy(_infoLine, "# User         : unknown \n");
#endif
}

// code by Max Moroz

const unsigned kMegaByte = 1024 * 1024;
const int kSmallChunks = 1000;
const unsigned kMaxAllocs = 20000;
const double MemUsageEps = 3;

// everything in bytes

inline long memused() {
  // cout << " Peak memory " << getPeakMemoryUsage() << endl;
  return static_cast<long>(getPeakMemoryUsage() * kMegaByte);
}

double getMemoryUsageEstimate() {
#if defined(_MSC_VER)
  return winMemUsage();
#endif

#if !defined(linux)
  return -1;
#endif
  static long prevMem = 0;
  static long extra;
  static int fail = 0;

  if (fail) return -1.;

  //      new_handler oldHndl;
  //      oldHndl=set_new_handler(mH);

  void *ptr[kMaxAllocs];
  unsigned numAllocs = 0;

  long last = memused();
  //       cout << "Memused : " << last << endl;
  if (last <= 0) return -1.;
  //      cerr<<"memused()="<<memused()<<endl;
  unsigned chunk = 8192;  // system allocates 8K chunks
  unsigned allocated = 0;
  int countSmallChunks = kSmallChunks;
  // if we allocate <8K, we'd get memused()+=8K, and allocated<8K - error
  while (1) {
    //              abkfatal(numAllocs<kMaxAllocs, "too many allocs");
    if (numAllocs >= kMaxAllocs) {
      abkwarn(0, "too many allocs (may be okay for 64-bit builds)");
      // FIFO destruction
      for (unsigned i = 0; i != numAllocs; ++i) free(ptr[i]);
      return -1.;
    }
    //              cerr<<"old: "<<memused()<<"; allocating "<<chunk<<"; now: ";
    if (!(ptr[numAllocs++] = malloc(chunk))) {
      fail = 1;
      return -1.;
    }
    //              cerr<<memused()<<endl;
    allocated += chunk;
    if (memused() > last + MemUsageEps) break;
    if (chunk <= kMegaByte && !countSmallChunks--) {
      chunk *= 2;
      countSmallChunks = kSmallChunks;
    }
  }

  //     int result=memused()-allocated-prevMem;

  /* LIFO destruction
          while (numAllocs) free(ptr[--numAllocs]); */

  // FIFO destruction
  for (unsigned i = 0; i != numAllocs; ++i) free(ptr[i]);

  extra = memused() - last;
  // handle extra correctly:
  // in some cases we need to add its prev value to current,
  // in some just store the new value

  prevMem = memused() - allocated;
  //      set_new_handler(oldHndl);
  return static_cast<double>(prevMem) / static_cast<double>(kMegaByte);
}

UserHomeDir::UserHomeDir() {
#if defined(linux)
  struct passwd *pwrec;
  pwrec = getpwuid(getuid());
  if (pwrec == NULL) {
    _infoLine = new char[8];
    strcpy(_infoLine, "unknown");
  } else {
    _infoLine = new char[strlen(pwrec->pw_dir) + 1];
    strcpy(_infoLine, pwrec->pw_dir);
  }
#elif defined(_MSC_VER)
  char *homedr = getenv("HOMEDRIVE");
  char *homepath = getenv("HOMEPATH");
  char *pathbuf = new char[_MAX_PATH];
  if (homedr == NULL || homepath == NULL)
    strcpy(pathbuf, "c:\\users\\default");
  else
    _makepath(pathbuf, homedr, homepath, NULL, NULL);
  unsigned len = strlen(pathbuf);
  if (pathbuf[len - 1] == '\\') pathbuf[len - 1] = '\0';
  _infoLine = new char[len + 1];  // +1 for the NUL strcpy() also copies
  strcpy(_infoLine, pathbuf);
  delete[] pathbuf;
#endif
}

ExecLocation::ExecLocation() {
  char buf[1000] = "";
#ifdef linux
  ssize_t pathLen = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
  if (pathLen > 0)
    buf[pathLen] = '\0';
  else
    buf[0] = '\0';
  int pos = strlen(buf) - 1;
  for (; pos >= 0; pos--) {
    if (buf[pos] == '/') {
      buf[pos] = '\0';
      break;
    }
  }
#endif

#if defined(_MSC_VER)
  ::GetModuleFileName(NULL, buf, 995);
  char drv[_MAX_DRIVE], dir[_MAX_DIR], fname[_MAX_FNAME], ext[_MAX_EXT];
  _splitpath(buf, drv, dir, fname, ext);
  _makepath(buf, drv, dir, "", "");
  int len = strlen(buf);
  if (buf[len - 1] == '\\') buf[len - 1] = '\0';
#endif

  _infoLine = new char[strlen(buf) + 1];
  strcpy(_infoLine, buf);
}

// <aaronnn>
static double getTotalPhysMemSize() {
#ifdef linux
#define AARON_BUFSZ 500
  char buf[AARON_BUFSZ];
  long totalPhysMem = -1;

  FILE *f;
  if ((f = fopen("/proc/meminfo", "r")) == NULL) return -1.0;
  if (fgets(buf, AARON_BUFSZ, f) == NULL) {
    fclose(f);
    return -1.0;
  }
  while (!feof(f) && !ferror(f)) {
    if (parseProcKilobytesLine(buf, "MemTotal:", &totalPhysMem)) break;

    // flush the pipe if buf was too small
    while (!feof(f) && !ferror(f) && strlen(buf) >= AARON_BUFSZ - 1)
      if (fgets(buf, AARON_BUFSZ, f) == NULL) break;

    if (fgets(buf, AARON_BUFSZ, f) == NULL) break;
  }
  fclose(f);

  if (totalPhysMem < 0) return -1.0;
  return double(totalPhysMem) / 1024.;
#undef AARON_BUFSZ
#else
  return -1.0;
#endif
}

static double getAvailablePhysMemSize() {
#ifdef linux
  struct stat s;
  if (stat("/proc/kcore", &s)) return -1.0;

  return s.st_size / (1024. * 1024.);
#else
  return -1.0;
#endif
}

static double getVmRSS() {
#ifdef linux
#define AARON_BUFSZ 500
  char buf[AARON_BUFSZ];
  long vmRSS = -1;

  FILE *f;
  if ((f = fopen("/proc/self/status", "r")) == NULL) return -1.0;
  if (fgets(buf, AARON_BUFSZ, f) == NULL) {
    fclose(f);
    return -1.0;
  }
  while (!feof(f) && !ferror(f)) {
    if (parseProcKilobytesLine(buf, "VmRSS:", &vmRSS)) break;

    // flush the pipe if buf was too small
    while (!feof(f) && !ferror(f) && strlen(buf) == AARON_BUFSZ - 1)
      if (fgets(buf, AARON_BUFSZ, f) == NULL) break;

    if (fgets(buf, AARON_BUFSZ, f) == NULL) break;
  }
  fclose(f);

  if (vmRSS < 0) return -1.0;
  return double(vmRSS) / 1024.;
#undef AARON_BUFSZ
#else
  return -1.0;
#endif
}

static double getMinorPageFaults() {
#ifdef linux
  struct rusage ru;

  if (getrusage(RUSAGE_SELF, &ru)) return -1.0;

  return ru.ru_minflt;
#else
  return -1.0;
#endif
}

static double getMajorPageFaults() {
#ifdef linux
  struct rusage ru;

  if (getrusage(RUSAGE_SELF, &ru)) return -1.0;

  return ru.ru_majflt;
#else
  return -1.0;
#endif
}

#ifdef linux
static void getProcSymlink(const char *elm, char *out, int outsz)
// find out where elm points to in /proc
// 'out' must have enough space allocated
{
  char path[MAXPATHLEN];

  if (snprintf(path, MAXPATHLEN, "/proc/self/%s", elm) >= MAXPATHLEN) {
    out[0] = '\0';
    return;
  }

  int written = readlink(path, out, outsz);

  if (written == -1)
    out[0] = '\0';
  else if (written == outsz)
    out[outsz - 1] = '\0';
  else
    out[written] = '\0';
}
#endif

#ifdef linux
static void getCPUInfo(char *out, int outsz)
/*
 * get cpu info and copy that in to buffer
 */
{
#define AARON_BUFSZ 500
  char buf[AARON_BUFSZ];
  char cpuID[256];
  double clock;
  unsigned numCPU = 1;
  char cpuInfo[AARON_BUFSZ];

  out[0] = '\0';

  FILE *f;
  if ((f = fopen("/proc/cpuinfo", "r")) == NULL) {
    out[0] = '\0';
    return;
  }
  bool modelNameFound, clockFound;
  modelNameFound = clockFound = false;
  if (fgets(buf, AARON_BUFSZ, f) == NULL) {
    fclose(f);
    return;
  }
  while (!feof(f) && !ferror(f)) {
    if (sscanf(buf, "model name %*s %255[-0-9_a-zA-z():$#+=. ]", cpuID) ==
        1) {  // (sizeof(cpuID) = 255+1)
      modelNameFound = true;
      clockFound = false;
    } else if (sscanf(buf, "cpu MHz %*s %lf", &clock) == 1) {
      clockFound = true;
    }

    if (modelNameFound && clockFound) {
      cpuID[255] = '\0';
      if (numCPU > 1)
        snprintf(cpuInfo, AARON_BUFSZ, "\n# CPU%d : %.2fMHz %s", numCPU, clock,
                 cpuID);
      else
        snprintf(cpuInfo, AARON_BUFSZ, "# CPU%d : %.2fMHz %s", numCPU, clock,
                 cpuID);
      cpuInfo[AARON_BUFSZ - 1] = '\0';
      strncat(out, cpuInfo, outsz - 1);
      out[outsz - 1] = '\0';
      modelNameFound = clockFound = false;
      numCPU++;
    }

    // flush the pipe if buf was too small
    while (!feof(f) && !ferror(f) && strlen(buf) == AARON_BUFSZ - 1)
      if (fgets(buf, AARON_BUFSZ, f) == NULL) break;

    if (fgets(buf, AARON_BUFSZ, f) == NULL) break;
  }
  fclose(f);
#undef AARON_BUFSZ
}
#endif

ExecInfo::ExecInfo()
//: Prints starting execution info
{
#ifdef linux
  std::string infoLine;
  char buf[MAXPATHLEN];

  infoLine += "# Working directory : ";
  getProcSymlink("cwd", buf, MAXPATHLEN);
  infoLine += buf;
  infoLine += "\n";

  infoLine += "# Executable : ";
  getProcSymlink(
      "exe", buf,
      MAXPATHLEN);  // WARNING: duplicate implementation in ExecLocation()
  infoLine += buf;
  infoLine += "\n";

  if (!gethostname(buf, MAXPATHLEN)) {
    infoLine += "# Hostname : ";
    infoLine += buf;
    infoLine += "\n";
  }

  getCPUInfo(buf, MAXPATHLEN);
  infoLine += buf;
  infoLine += "\n";

  _infoLine = new char[infoLine.size() + 1];
  strcpy(_infoLine, infoLine.c_str());
#else
  _infoLine = new char[1];
  _infoLine[0] = '\0';
#endif
}
// </aaronnn>
