
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __RealSimTk_H_
#define __RealSimTk_H__

#ifndef RealType
#define RealType 1
#endif 

#if RealType == 1 

#define RealOpenMM float
#define SQRT       sqrtf
#define POW        powf
#define SIN        sinf
#define COS        cosf
#define TAN        tanf

// LOG is used in Vishal's gpu code; modifying LOG -> LN 
#define LN         logf

#define EXP        expf
#define FABS       fabsf
#define ACOS       acosf
#define ASIN       asinf
#define ATAN       atanf
#define TANH       tanhf

#define ATOF       atoff

#define PI_M             3.141592653589f
#define TWO_SIX          1.122462048309372981f
#define RADIAN           57.29577951308f
#define LOG_TEN          2.302585092994045684f
#define SQRT_TWO         1.41421356237309504f
#define RADIAN_INVERSE   0.01745329252f

#else

#define RealOpenMM double
#define SQRT       sqrt
#define POW        pow
#define SIN        sin
#define COS        cos
#define TAN        tan

// LOG is used in Vishal's gpu code; modifying LOG -> LN 
#define LN         log

#define EXP        exp
#define FABS       fabs
#define ACOS       acos
#define ASIN       asin
#define ATAN       atan
#define TANH       tanh

#define ATOF       atof

#define PI              3.141592653589
#define TWO_SIX         1.122462048309372981
#define RADIAN         57.29577951308
#define LOG_TEN         2.302585092994045684
#define SQRT_TWO        1.41421356237309504
#define RADIAN_INVERSE  0.01745329252

#endif

#define DOT3(u,v) ((u[0])*(v[0]) + (u[1])*(v[1]) + (u[2])*(v[2]))

#define MATRIXDOT3(u,v) u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + \
                        u[3]*v[3] + u[4]*v[4] + u[5]*v[5] + \
                        u[6]*v[6] + u[7]*v[7] + u[8]*v[8]


#endif // __RealSimTk_H__
