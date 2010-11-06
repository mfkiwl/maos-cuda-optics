/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __DRAWDAEMON_H
#define __DRAWDAEMON_H
#include "misc.h"
#define CSIZE       512
enum{
    FIFO_START,
    FIFO_DATA,
    FIFO_SHM,
    FIFO_POINTS,
    FIFO_STYLE,
    FIFO_CIRCLE,
    FIFO_LIMIT,
    FIFO_FIG,
    FIFO_NAME,
    FIFO_TITLE,
    FIFO_XLABEL,
    FIFO_YLABEL,
    FIFO_MAXMIN,
    FIFO_END
};

#ifndef CAIRO_FORMAT_A8
#define CAIRO_FORMAT_RGB24 0x01
#define CAIRO_FORMAT_A8 0x02
#endif
#define T_CHAR   (0x1<<16|CAIRO_FORMAT_A8<<8    | sizeof(char))
#define T_INT    (0x2<<16|CAIRO_FORMAT_RGB24<<8 | sizeof(int))

#define T_DOUBLE (0x3<<16 | sizeof(double))
#define T_LONG   (0x4<<16 | sizeof(long))
#define T_DCOMPLEX (0x5<<16 | sizeof(dcomplex))
#define DEFAULT_FONT_SIZE 12
#endif