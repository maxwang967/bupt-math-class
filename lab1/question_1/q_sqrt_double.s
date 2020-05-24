# 体系结构 x86_64
my_sqrt(double):
        movapd  xmm2, xmm0
        pxor    xmm0, xmm0
        comisd  xmm0, xmm2
        ja      .L4
        movsd   xmm5, QWORD PTR .LC2[rip]
        movq    xmm4, QWORD PTR .LC3[rip]
        movapd  xmm0, xmm2
        movsd   xmm3, QWORD PTR .LC4[rip]
.L3:
        movapd  xmm1, xmm2
        divsd   xmm1, xmm0
        addsd   xmm0, xmm1
        mulsd   xmm0, xmm5
        movapd  xmm1, xmm0
        mulsd   xmm1, xmm0
        subsd   xmm1, xmm2
        andpd   xmm1, xmm4
        comisd  xmm1, xmm3
        ja      .L3
        ret
.L4:
        movsd   xmm0, QWORD PTR .LC0[rip]
        ret
.LC0:
        .long   0
        .long   -1074790400
.LC2:
        .long   0
        .long   1071644672
.LC3:
        .long   -1
        .long   2147483647
        .long   0
        .long   0
.LC4:
        .long   -1629006314
        .long   1020396463