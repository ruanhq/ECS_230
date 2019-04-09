/*
 * timing2.c
 *
 * test readTSC function
 *
 * compile using: gcc timing2.c -lm
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

long long readTSC(void)
{
 /* read the time stamp counter on Intel x86 chips */
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
            : "=mr" (ticks.part[0]),
              "=mr" (ticks.part[1])
            : /* no inputs */
            : "eax", "edx");
  return ticks.complete;
}

int main(int argc, char** argv)
{
  int n = atoi(argv[1]);
  long long tm;
  double x;
  tm = readTSC();
  x = sqrt(10.0);
  tm = readTSC() - tm;
  printf("%lld TSC clocks\n",tm);
  printf("%f\n",x);


  double sum = 0.0;
  int i;
  tm = readTSC();
  for ( i = 0; i < n; i++ )
    sum += sqrt((double) i);
  tm = readTSC() - tm;
  printf("%lld TSC clocks\n",tm);
  printf("%15.6e\n",sum);

  return 0;
}

