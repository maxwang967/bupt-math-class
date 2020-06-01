# 体系结构 intel x86_64 
my_sqrt:
        vxorps  xmm1, xmm1, xmm1        # if (c < 0)
        vcomiss xmm1, xmm0
        ja      .L5
        vmovd   eax, xmm0               # int exp = (ix >> 23) & 0xFF;
        mov     ecx, 23
        vxorps  xmm1, xmm1, xmm1
        vcvtss2sd       xmm0, xmm0, xmm0
        sarx    eax, eax, ecx
        movzx   eax, al
        sub     eax, 127                # int log2 = exp - 127;
        mov     edx, eax                # double h = (long long)1 << (e / 2);
        shr     edx, 31
        add     eax, edx
        mov     edx, 1
        sar     eax
        shlx    rax, rdx, rax
        vcvtsi2sd       xmm1, xmm1, rax
        vmovsd  xmm2, xmm1, xmm1
        vmulsd  xmm1, xmm0, QWORD PTR p[rip]        # t = p * c / h + q * h;
        vdivsd  xmm1, xmm1, xmm2    
        vfmadd231sd     xmm1, xmm2, QWORD PTR q[rip]
        vdivsd  xmm2, xmm0, xmm1         # t = (c / t + t) / 2.0;
        vaddsd  xmm1, xmm2, xmm1
        vmovsd  xmm2, QWORD PTR .LC2[rip]
        vmulsd  xmm1, xmm1, xmm2 
        vdivsd  xmm0, xmm0, xmm1        # t = (c / t + t) / 2.0;
        vaddsd  xmm0, xmm0, xmm1
        vmulsd  xmm0, xmm0, xmm2 
        vcvtsd2ss       xmm0, xmm0, xmm0 
        ret     
.L5:
        vmovss  xmm0, DWORD PTR .LC0[rip]
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
        .long   -1082130432
.LC2:
        .long   0
        .long   1071644672