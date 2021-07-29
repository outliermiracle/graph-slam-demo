/*
graph-slam-demo
a simple SLAM example demonstrating the usage of graph optimization and Lie algebra
by outliermiracle

compile: g++ -o loop_closing loop_closing.cpp
pre-request: eigen v3.2.92
usage: ./loop_closing
output: file data.csv, with optimized points 
*/

#include <iostream>
#include <Eigen/Dense>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace Eigen;

Vector4f x[5];
Matrix4f T[5];
Vector4f e[5];

const Matrix3f I3 = Matrix3f::Identity();

Vector3f t[5];
Matrix3f R[5];

Matrix<float, 6, 1> ksi[5];
Matrix<float, 3, 6> LR[5];

Matrix3f make_antisymmetric(Vector3f a)
{
    Matrix3f a_antisymmetric = Matrix3f::Zero();
    a_antisymmetric(1,0) = a(2);
    a_antisymmetric(0,1) = -a(2);
    a_antisymmetric(2,0) = -a(1);
    a_antisymmetric(0,2) = a(1);
    a_antisymmetric(2,1) = a(0);
    a_antisymmetric(1,2) = -a(0);
    return a_antisymmetric;
}

int main()
{
    FILE *fp;
    fp = fopen("data.csv", "w");
    int i, j;

    t[1] << 0, 1, 0;
    t[2] << 1+0.1, 0+0.1, 0;
    t[3] << 0, -1, 0;
    t[4] << -1, 0, 0;
    t[0] << 0, 0, 0;

    for(i=0; i<5; i++)
    {
        T[i].block(0,0,3,3) = I3;
        T[i].block(0,3,3,1) = t[i];
        T[i].block(3,0,1,4) << 0,0,0,1;
        printf("T[%d] = \n", i);
        cout << T[i] << endl << endl;
    }

    x[0].segment(0, 3) = t[0];
    x[0].segment(3, 1) << 1;
    for(i=1; i<5; i++)
    {
        x[i] = T[i] * x[i-1];
    }

    for(i=0; i<5; i++)
    {
        printf("x[%d] = ", i);
        cout << x[i].transpose().block(0,0,1,3) << endl;
        fprintf(fp, "%f, %f\n", x[i](0), x[i](1));
    }
    fprintf(fp, "\n");

    //main loop
    for(int loop_idx=0; loop_idx<100; loop_idx++)
    {
        printf("==========loop %d==========\n", loop_idx);
        Matrix<float, 6, 1> b;
        b.setZero();
        Matrix<float, 6, 6> H;
        H.setZero();
        float F = 0;
        for(i=0; i<5; i++)
        {
            e[i] = T[(i+1)%5] * x[i] - x[(i+1)%5];
            F += e[i].squaredNorm();
            R[i] = T[i].block(0,0,3,3);
            float tr_R = R[i].trace();
            float theta = acos((tr_R - 1) / 2);

            EigenSolver<Matrix3f> eig(R[i]);
            Array<std::complex<float>, Dynamic, 1> eigen_values = eig.eigenvalues();

            int size = eigen_values.size();
            bool found_eigen = false;
            int index_eigen;
            for(j=0; j<size; j++)
            {
                if(eigen_values[j].imag() == 0 && eigen_values[j].real() == 1)
                {
                    index_eigen = j;
                    found_eigen = true;
                    break;
                }
            }

            if(!found_eigen)
            {
                printf("eigen value error!\n");
                printf("R[%d].eigenvalues() = \n", i);
                cout << R[i].eigenvalues() << endl << endl;
                EigenSolver<Matrix3f> eig(R[i]);
                printf("R[%d].eigenvectors() = \n", i);
                cout << eig.eigenvectors() << endl << endl;
                return -1;
            }

            Vector3f a = eig.eigenvectors().real().col(index_eigen);
            Vector3f phi = theta * a;
            Matrix3f aaT = a * a.transpose();
            Matrix3f a_antisymmetric = make_antisymmetric(a);

            float s_theta = sin(theta);
            float c_theta = cos(theta);
            float f_theta;
            if(fabs(theta) > 0.001) f_theta = s_theta / theta;
            else f_theta = 1;
            Matrix3f J = f_theta * I3;
            J += (1 - f_theta) * aaT;
            if(fabs(theta) > 0.001) J += (1 - c_theta) / theta * a_antisymmetric;

            //线性方程求解 J * rho = t[i];
            Vector3f rho = J.colPivHouseholderQr().solve(t[i]);

            ksi[i].block(0,0,3,1) = phi;
            ksi[i].block(3,0,3,1) = rho;

            Matrix3f x_antisymmetric = make_antisymmetric(x[i].block(0,0,3,1));

            Matrix<float, 3, 6> Lk;
            Lk.block(0,0,3,3) = -x_antisymmetric;
            Lk.block(0,3,3,3) = I3;
            LR[i] = Lk;
            Lk = R[(i+1)%5] * Lk;

            Matrix<float, 1, 3> ekT = e[i].transpose().block(0,0,1,3);

            Matrix<float, 6, 1> bk;
            bk = (ekT * Lk).transpose();
            b += bk;
            Matrix<float, 6, 6> Hk;
            Hk = Lk.transpose() * Lk;
            H += Hk;
        }

        //线性方程求解 H * d_ksi = -b;
        Matrix<float, 6, 1> d_ksi;
        d_ksi = H.colPivHouseholderQr().solve(-b) * 0.1;

        for(i=1; i<5; i++)
        {
            ksi[i] += d_ksi;
            Vector3f phi = ksi[i].block(0,0,3,1);
            Vector3f rho = ksi[i].block(3,0,3,1);
            float theta = phi.norm();
            Vector3f a;
            if(fabs(theta) > 0.01) a = phi / theta;
            else a << 1,0,0;
            Matrix3f aaT = a * a.transpose();
            Matrix3f a_antisymmetric = make_antisymmetric(a);

            float s_theta = sin(theta);
            float c_theta = cos(theta);
            float f_theta;
            if(fabs(theta) > 0.01) f_theta = s_theta / theta;
            else f_theta = 1;
            Matrix3f J = f_theta * I3;
            J += (1 - f_theta) * aaT;
            if(fabs(theta) > 0.01) J += (1 - c_theta) / theta * a_antisymmetric;

            R[i] = c_theta * I3 + (1 - c_theta) * aaT + s_theta * a_antisymmetric;
            T[i].block(0,0,3,3) = R[i];
            t[i] = J * rho;
            T[i].block(0,3,3,1) = t[i];
            T[i].block(3,0,1,4) << 0,0,0,1;
        }
        for(i=1; i<5; i++)
        {
            x[i] = T[i] * x[i-1];
        }
        for(i=0; i<5; i++)
        {
            printf("x[%d] = ", i);
            cout << x[i].transpose().block(0,0,1,3) << endl;
            if(loop_idx % 20 == 0) fprintf(fp, "%f, %f\n", x[i](0), x[i](1));
        }
        if(loop_idx % 20 == 0) fprintf(fp, "\n");
        printf("error = %f\n", F);
    }
    fclose(fp);

    return 0;
}
