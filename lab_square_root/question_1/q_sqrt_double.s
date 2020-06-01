# 体系结构 intel x86_64 
my_sqrt_d:
        vxorpd  xmm1, xmm1, xmm1        # if (c < 0)
        vcomisd xmm1, xmm0
        ja      .L5
        vmovq   rax, xmm0               # long long exp = (ix >> 52) & 0xFFF;
        mov     ecx, 52
        vxorps  xmm1, xmm1, xmm1
        shrx    rax, rax, rcx
        lea     rdx, [rax-1023]         # long long log2 = exp - 1023;
        mov     eax, edx                # double h = (long long)1 << (e / 2);
        shr     eax, 31
        add     eax, edx
        mov     edx, 1
        sar     eax
        shlx    rax, rdx, rax
        vcvtsi2sd       xmm1, xmm1, rax
        vmovsd  xmm2, xmm1, xmm1
        vmulsd  xmm1, xmm0, QWORD PTR p[rip]    # double t = p * c / h + q * h;
        vdivsd  xmm1, xmm1, xmm2
        vfmadd231sd     xmm1, xmm2, QWORD PTR q[rip]
        vdivsd  xmm2, xmm0, xmm1        # t = (c / t + t) / 2.0;
        vaddsd  xmm1, xmm2, xmm1
        vmulsd  xmm2, xmm1, xmm3
        vdivsd  xmm1, xmm0, xmm2        # t = (c / t + t) / 2.0;
        vaddsd  xmm1, xmm1, xmm2
        vmulsd  xmm1, xmm1, xmm3
        vdivsd  xmm0, xmm0, xmm1        # t = (c / t + t) / 2.0;
        vaddsd  xmm0, xmm0, xmm1
        vmulsd  xmm0, xmm0, xmm3
        ret
.L5:
        vmovsd  xmm0, QWORD PTR .LC0[rip]   # <retval>,
        ret     
sqrt_2:
        .long   1719614413
        .long   1073127582
q:
        .long   857772838
        .long   1071843209
p:
        .long   -1711476940
        .long   1071284857
.LC0:
        .long   0
        .long   -1074790400
.LC2:
        .long   0
        .long   1071644672