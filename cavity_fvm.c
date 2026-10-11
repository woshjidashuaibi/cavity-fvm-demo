/* cavity_simple.c
 * 2D Lid-Driven Cavity, FVM + SIMPLE, Colocated grid, pure C, no libs
 * Re = 100
 * Compile: gcc cavity_simple.c -o cavity -lm
 * Run: ./cavity
 * Output: cavity.dat (Tecplot format)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NI 40          // x方向网格数（单元）
#define NJ 40          // y方向网格数（单元）
#define MAX_ITER 8000  // 最大迭代步数
#define TOL 1e-6       // 连续性残差收敛阈值
#define RE 100.0       // Reynolds number
#define U_TOP 1.0      // 顶盖速度

double dx = 1.0 / NI;
double dy = 1.0 / NJ;
double nu = U_TOP * 1.0 / RE;  // kinematic viscosity

// 变量：同位网格，单元中心存储
double u[NI+2][NJ+2];
double v[NI+2][NJ+2];
double p[NI+2][NJ+2];
double pprime[NI+2][NJ+2];
double ap[NI+2][NJ+2], ae[NI+2][NJ+2], aw[NI+2][NJ+2], an[NI+2][NJ+2], as[NI+2][NJ+2];
double b[NI+2][NJ+2];

// 松弛因子
const double alpha_u = 0.7;
const double alpha_p = 0.3;

// 求解五对角方程组（点迭代，Jacobi，简单易实现，适合教学）
void solve_uv(double var[NI+2][NJ+2])
{
    int it, i, j;
    double old;
    for(it=0;it<150;it++)
    {
        for(i=1;i<=NI;i++)
        for(j=1;j<=NJ;j++)
        {
            old = var[i][j];
            var[i][j] = ( ae[i][j]*var[i+1][j] + aw[i][j]*var[i-1][j]
                        + an[i][j]*var[i][j+1] + as[i][j]*var[i][j-1] + b[i][j] ) / ap[i][j];
            var[i][j] = old + alpha_u*(var[i][j] - old);
        }
    }
}

// 压力修正方程求解
void solve_pp()
{
    int it, i, j;
    double old;
    for(it=0;it<200;it++)
    {
        for(i=1;i<=NI;i++)
        for(j=1;j<=NJ;j++)
        {
            old = pprime[i][j];
            pprime[i][j] = ( ae[i][j]*pprime[i+1][j] + aw[i][j]*pprime[i-1][j]
                            + an[i][j]*pprime[i][j+1] + as[i][j]*pprime[i][j-1] + b[i][j] ) / ap[i][j];
            pprime[i][j] = old + alpha_p*(pprime[i][j] - old);
        }
    }
}

// 组装u动量方程系数
void assemble_u()
{
    int i,j;
    double Fw,Fe,Fs,Fn;
    double D = nu/dx;
    for(i=1;i<=NI;i++)
    for(j=1;j<=NJ;j++)
    {
        // 界面对流通量（一阶迎风格式）
        Fw = 0.5*(u[i-1][j]+u[i][j]) * dy;
        Fe = 0.5*(u[i][j]+u[i+1][j]) * dy;
        Fs = 0.5*(v[i][j-1]+v[i][j]) * dx;
        Fn = 0.5*(v[i][j]+v[i][j+1]) * dx;

        aw[i][j] = D + fmax(Fw, 0.0);
        ae[i][j] = D + fmax(-Fe,0.0);
        as[i][j] = nu/dy + fmax(Fs,0.0);
        an[i][j] = nu/dy + fmax(-Fn,0.0);
        ap[i][j] = aw[i][j]+ae[i][j]+as[i][j]+an[i][j];
        b[i][j] = (p[i-1][j] - p[i][j]) * dy;
    }
    // u边界条件
    for(j=0;j<=NJ+1;j++) u[0][j]=0.0, u[NI+1][j]=0.0;
    for(i=0;i<=NI+1;i++) u[i][0]=0.0, u[i][NJ+1]=U_TOP;
}

// 组装v动量方程系数
void assemble_v()
{
    int i,j;
    double Fw,Fe,Fs,Fn;
    double D = nu/dy;
    for(i=1;i<=NI;i++)
    for(j=1;j<=NJ;j++)
    {
        Fw = 0.5*(u[i-1][j]+u[i][j]) * dy;
        Fe = 0.5*(u[i][j]+u[i+1][j]) * dy;
        Fs = 0.5*(v[i][j-1]+v[i][j]) * dx;
        Fn = 0.5*(v[i][j]+v[i][j+1]) * dx;

        aw[i][j] = nu/dx + fmax(Fw, 0.0);
        ae[i][j] = nu/dx + fmax(-Fe,0.0);
        as[i][j] = D + fmax(Fs, 0.0);
        an[i][j] = D + fmax(-Fn,0.0);
        ap[i][j] = aw[i][j]+ae[i][j]+as[i][j]+an[i][j];
        b[i][j] = (p[i][j-1] - p[i][j]) * dx;
    }
    // v边界条件，所有壁面v=0
    for(j=0;j<=NJ+1;j++) v[0][j]=0.0, v[NI+1][j]=0.0;
    for(i=0;i<=NI+1;i++) v[i][0]=0.0, v[i][NJ+1]=0.0;
}

// 组装压力修正方程
double assemble_pp()
{
    int i,j;
    double res=0.0;
    double du, dv;
    for(i=1;i<=NI;i++)
    for(j=1;j<=NJ;j++)
    {
        double d_u = dy / ap[i][j];
        double d_v = dx / ap[i][j];
        aw[i][j] = d_u;
        ae[i][j] = d_u;
        as[i][j] = d_v;
        an[i][j] = d_v;
        ap[i][j] = aw[i][j]+ae[i][j]+as[i][j]+an[i][j];
        // 连续性源项 b = mass imbalance
        du = dy*(u[i][j]-u[i+1][j]);
        dv = dx*(v[i][j]-v[i][j+1]);
        b[i][j] = du + dv;
        res += fabs(b[i][j]);
    }
    // 压力边界：p'在一个点固定为0（参考压力）
    pprime[1][1]=0.0;
    return res/(NI*NJ);
}

// 速度修正
void correct_uv()
{
    int i,j;
    for(i=1;i<=NI;i++)
    for(j=1;j<=NJ;j++)
    {
        double du = dy/ap[i][j];
        double dv = dx/ap[i][j];
        u[i][j] -= du*(pprime[i+1][j]-pprime[i][j]);
        v[i][j] -= dv*(pprime[i][j+1]-pprime[i][j]);
    }
    // 修正压力
    for(i=1;i<=NI;i++)
    for(j=1;j<=NJ;j++)
    {
        p[i][j] += alpha_p * pprime[i][j];
        pprime[i][j] = 0.0;
    }
    // BC再次强制
    for(j=0;j<=NJ+1;j++) u[0][j]=0.0, u[NI+1][j]=0.0;
    for(i=0;i<=NI+1;i++) u[i][0]=0.0, u[i][NJ+1]=U_TOP;
    for(j=0;j<=NJ+1;j++) v[0][j]=0.0, v[NI+1][j]=0.0;
    for(i=0;i<=NI+1;i++) v[i][0]=0.0, v[i][NJ+1]=0.0;
}

// 输出Tecplot格式
void write_tecplot()
{
    FILE *fp = fopen("cavity.dat","w");
    int i,j;
    fprintf(fp,"TITLE=\"Lid Driven Cavity Re=100\"\n");
    fprintf(fp,"VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"P\"\n");
    fprintf(fp,"ZONE I=%d, J=%d, F=POINT\n",NI,NJ);
    for(j=1;j<=NJ;j++)
    for(i=1;i<=NI;i++)
    {
        double x = (i-0.5)*dx;
        double y = (j-0.5)*dy;
        fprintf(fp,"%lf %lf %lf %lf %lf\n",x,y,u[i][j],v[i][j],p[i][j]);
    }
    fclose(fp);
    printf("Result saved to cavity.dat\n");
}

int main()
{
    int iter;
    double res_cont;
    // init zero
    for(int i=0;i<=NI+1;i++)
    for(int j=0;j<=NJ+1;j++)
    {
        u[i][j]=0;v[i][j]=0;p[i][j]=0;pprime[i][j]=0;
    }

    printf("SIMPLE start, Re=%g, mesh %dx%d\n",RE,NI,NJ);
    for(iter=0;iter<MAX_ITER;iter++)
    {
        // 1. Solve momentum for u*
        assemble_u();
        solve_uv(u);
        // 2. Solve momentum for v*
        assemble_v();
        solve_uv(v);
        //3. Pressure correction
        res_cont = assemble_pp();
        solve_pp();
        //4. Correct u,v,p
        correct_uv();

        if(iter%100 == 0)
        {
            printf("Iter %d, Continuity residual = %e\n",iter,res_cont);
        }
        if(res_cont < TOL)
        {
            printf("Converged at iter=%d, residual=%e\n",iter,res_cont);
            break;
        }
    }
    write_tecplot();
    return 0;
}
