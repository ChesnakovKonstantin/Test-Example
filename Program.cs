using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Drawing;

namespace ImageReadCS
{
    class Program
    {
        static ColorFloatImage Gauss(ColorFloatImage image, int n, double sigma)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width, image.Height), tmp = new ColorFloatImage(image.Width + 2 * n, image.Height + 2 * n);
            float[,] masG = new float[(2 * n + 1), (2 * n + 1)];
            //double n1 = 0, m1 = 0, s = 0;
            //double res = 0;
            //bool t = true;
            //Random rnd = new Random();
            //float[] masr = new float[(2 * n + 1) * (2 * n + 1)], masb = new float[(2 * n + 1) * (2 * n + 1)], masg = new float[(2 * n + 1) * (2 * n + 1)];
            int k1 = 2 * n + 1;
            float coeff = (float)(1.0f / (2.0f * Math.PI * sigma * sigma));
            float sum = 0;
            //float[,] A = new float[k, k];
            for (int i = 0; i < 2 * n + 1; i++)
                for (int j = 0; j < 2 * n + 1; j++)
                {
                    masG[i, j] = coeff * (float)Math.Exp(-(Math.Pow(n - i, 2) + Math.Pow(n - j, 2)) / (2 * sigma * sigma));
                    sum += masG[i, j];
                    Console.WriteLine(masG[i,j]);
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[n + j, n + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < n; j++)
                {
                    tmp[j, i+n] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < n; j++)
                {
                    tmp[j + n + image.Width, i+n] = image[image.Width - 1, i];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + n, i] = image[j, 0];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + n, i + n + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmp[j, i] = image[0, 0];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmp[j + n + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmp[j + n + image.Width, i + n + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmp[j, i + n + image.Height] = image[0, image.Height - 1];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    ColorFloatPixel p = image[j, i];
                    p.r = 0;
                    p.g = 0;
                    p.b = 0;
                    for (int k = 0; k < 2 * n + 1; k++)
                        for (int m = 0; m < 2 * n + 1; m++)
                        {
                            ColorFloatPixel p1 = tmp[m + j, k + i];
                            p.r = p.r + p1.r * masG[k, m]/sum;
                            p.g = p.g + p1.g * masG[k, m]/sum;
                            p.b = p.b + p1.b * masG[k, m]/sum;
                        }
                    image1[j, i] = p;
                }
            return image1;
        }

        static float BiLinearInterp(float y12, float y22, float y11, float y21, int x1, int y1, int n, int x, int y)
        {
            return /*(1/((n)^2))**/((-y11)*(x1+n-x)*(y1-n-y) + (-y21)*(x-x1)*(y1-n-y) + (-y12)*(x1+n-x)*(y-y1) + (-y22)*(x-x1)*(y-y1))/(n*n);
        }

        static ColorFloatImage BileanIntInt(ColorFloatImage image, int n)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width*n, image.Height*n),
                            tmp = new ColorFloatImage(image.Width* n + 1, image.Height*n + 1);
             /*for (int y = 0; y < image.Height; y++)
                 for (int x = 0; x < image.Width; x++)
                 {
                     image1[x * n, y * n] = image[x, y];
                 }*/
             for (int i = 0; i < image.Height; i++)
                 for (int j = 0; j < image.Width; j++)
                 {
                     tmp[j * n, i * n] = image[j, i];
                 }
             for (int i = 0; i < image.Height; i++)
                 {
                     tmp[image.Width*n, i*n] = image[image.Width - 1, i];
                 }
             for (int j = 0; j < image.Width; j++)
                 {
                     tmp[j*n, image.Height*n] = image[j, image.Height - 1];
                 }
             tmp[image.Width*n, image.Height*n] = image[image.Width - 1, image.Height - 1];
                 
            for (int y = 0; y < image.Height; y++)
                for (int x = 0; x < image.Width; x++)
                {
                    ColorFloatPixel p = image[x, y];
                    ColorFloatPixel p11 = tmp[x * n, y * n + n];
                    ColorFloatPixel p12 = tmp[x * n, y * n];
                    ColorFloatPixel p21 = tmp[x * n + n, y * n + n];
                    ColorFloatPixel p22 = tmp[x * n + n, y * n];
                    for (int j = 1; j < n; j++)
                        for (int i = 1; i < n; i++)
                        {
                            p.r = BiLinearInterp(p12.r, p22.r, p11.r, p21.r, x * n, y * n + n, n, x * n + i, y * n + j);
                            p.b = BiLinearInterp(p12.b, p22.b, p11.b, p21.b, x * n, y * n + n, n, x * n + i, y * n + j);
                            p.g = BiLinearInterp(p12.g, p22.g, p11.g, p21.g, x * n, y * n + n, n, x * n + i, y * n + j);
                            tmp[x * n + i, y * n + j] = p;
                        }
                    for (int i = 1; i < n; i++)
                    {
                        p.r = BiLinearInterp(p12.r, p22.r, p11.r, p21.r, x * n, y * n + n, n, x * n + i, y * n);
                        p.b = BiLinearInterp(p12.b, p22.b, p11.b, p21.b, x * n, y * n + n, n, x * n + i, y * n);
                        p.g = BiLinearInterp(p12.g, p22.g, p11.g, p21.g, x * n, y * n + n, n, x * n + i, y * n);
                        tmp[x * n + i, y * n] = p;
                    }
                    for (int i = 1; i < n; i++)
                    {
                        p.r = BiLinearInterp(p12.r, p22.r, p11.r, p21.r, x * n, y * n + n, n, x * n, y * n + i);
                        p.b = BiLinearInterp(p12.b, p22.b, p11.b, p21.b, x * n, y * n + n, n, x * n, y * n + i);
                        p.g = BiLinearInterp(p12.g, p22.g, p11.g, p21.g, x * n, y * n + n, n, x * n, y * n + i);
                        tmp[x * n, y * n + i] = p;
                    }
                }
            for (int y = 0; y < image.Height*n; y++)
                for (int x = 0; x < image.Width * n; x++)
                {
                    image1[x, y] = tmp[x, y];
                }
            return image1;
        }

        static float BiCubeInterp(float [] mas, int x1, int y1, int n, float x, float y)
        {
            return (mas[0] * (x - n) * (x - 2 * n) * (x + n) * (y - n) * (y - 2 * n) * (y + n) / (4) +
                   mas[1] * x * (x - 2 * n) * (x + n) * (y - n) * (y - 2 * n) * (y + n) / (-4) +
                   mas[2] * (x - n) * (x - 2 * n) * (x + n) * (y) * (y - 2 * n) * (y + n) / (-4) +
                   mas[3] * (x) * (y) * (x - 2 * n) * (x + n) * (y - 2 * n) * (y + n) / 4 +
                   mas[4] * (x - n) * (y - n) * (x - 2 * n) * (x) * (y - 2 * n) * (y + n) / (-12) +
                   mas[5] * (x - n) * (x - 2 * n) * (x + n) * y * (y - n) * (y - 2 * n) / (-12) +
                   mas[6] * (x - n) * (x - 2 * n) * x * y * (y - 2 * n) * (y + n) / 12 +
                   mas[7] * (x - 2 * n) * x * (x + n) * y * (y - n) * (y - 2 * n) / 12 +
                   mas[8] * (x - n) * x * (x + n) * (y - n) * (y - 2 * n) * (y + n) / 12 +
                   mas[9] * (x - n) * (x - 2 * n) * (x + n) * y * (y - n) * (y + n) / 12 +
                   mas[10] * (x - n) * (x - 2 * n) * x * y * (y - n) * (y - 2 * n) / 36 +
                   mas[11] * (x - n) * x * (x + n) * y * (y - 2 * n) * (y + n) / (-12) +
                   mas[12] * (x - 2 * n) * x * (x + n) * y * (y - n) * (y + n) / (-12) +
                   mas[13] * (x - n) * x * (x + n) * y * (y - n) * (y - 2 * n) / (-36) +
                   mas[14] * (x - n) * (x - 2 * n) * x * y * (y - n) * (y + n) / (-36) +
                   mas[15] * (x - n) * x * (x + n) * y * (y - n) * (y + n) / (36)) / (n * n * n * n * n * n);
        }

        static ColorFloatImage BicubeIntInt1(ColorFloatImage image, int n)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width * n, image.Height * n),
                            tmp = new ColorFloatImage(image.Width + 5, image.Height + 5);
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[1 + j, 1 + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j,1 + i] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 4; j++)
                {
                    tmp[j + 1 + image.Width, i+1] = image[image.Width - 1, i];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 1, i] = image[j, 0];
                }
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 1, i + 1 + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, i] = image[0, 0];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < 4; j++)
                {
                    tmp[j + 1 + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                {
                    tmp[j + 1 + image.Width, i + 1 + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, i + 1 + image.Height] = image[0, image.Height - 1];
                }
            float[] masr = new float[16];
            float[] masb = new float[16];
            float[] masg = new float[16];
            for (int y = 0; y < image.Height; y= y + 1)
                for (int x = 0; x < image.Width; x = x + 1)
                {
                    ColorFloatPixel p = image[0, 0];
                    ColorFloatPixel p00 = tmp[x + 1, y + 1];
                    masr[0] = p00.r;
                    masb[0] = p00.b;
                    masg[0] = p00.g;
                    ColorFloatPixel p01 = tmp[x + 1, y + 2];
                    masr[1] = p01.r;
                    masb[1] = p01.b;
                    masg[1] = p01.g;
                    ColorFloatPixel p10 = tmp[x + 2, y + 1];
                    masr[2] = p10.r;
                    masb[2] = p10.b;
                    masg[2] = p10.g;
                    ColorFloatPixel p11 = tmp[x + 2, y + 2];
                    masr[3] = p11.r;
                    masb[3] = p11.b;
                    masg[3] = p11.g;
                    ColorFloatPixel p0_1 = tmp[x + 1, y];
                    masr[4] = p0_1.r;
                    masb[4] = p0_1.b;
                    masg[4] = p0_1.g;
                    ColorFloatPixel p_10 = tmp[x, y + 1];
                    masr[5] = p_10.r;
                    masb[5] = p_10.b;
                    masg[5] = p_10.g;
                    ColorFloatPixel p1_1 = tmp[x + 2, y];
                    masr[6] = p1_1.r;
                    masb[6] = p1_1.b;
                    masg[6] = p1_1.g;
                    ColorFloatPixel p_11 = tmp[x, y + 2];
                    masr[7] = p_11.r;
                    masb[7] = p_11.b;
                    masg[7] = p_11.g;
                    ColorFloatPixel p02 = tmp[x + 1, y + 3];
                    masr[8] = p02.r;
                    masb[8] = p02.b;
                    masg[8] = p02.g;
                    ColorFloatPixel p20 = tmp[x + 3, y + 1];
                    masr[9] = p20.r;
                    masb[9] = p20.b;
                    masg[9] = p20.g;
                    ColorFloatPixel p_1_1 = tmp[x, y];
                    masr[10] = p_1_1.r;
                    masb[10] = p_1_1.b;
                    masg[10] = p_1_1.g;
                    ColorFloatPixel p12 = tmp[x + 2, y + 3];
                    masr[11] = p12.r;
                    masb[11] = p12.b;
                    masg[11] = p12.g;
                    ColorFloatPixel p21 = tmp[x + 3, y + 2];
                    masr[12] = p21.r;
                    masb[12] = p21.b;
                    masg[12] = p21.g;
                    ColorFloatPixel p_12 = tmp[x, y + 3];
                    masr[13] = p_12.r;
                    masb[13] = p_12.b;
                    masg[13] = p_12.g;
                    ColorFloatPixel p2_1 = tmp[x + 3, y];
                    masr[14] = p2_1.r;
                    masb[14] = p2_1.b;
                    masg[14] = p2_1.g;
                    ColorFloatPixel p22 = tmp[x + 3, y + 3];
                    masr[15] = p22.r;
                    masb[15] = p22.b;
                    masg[15] = p22.g;
                    if (y == 0 && x == 2)
                        for (int k = 0; k < 16; k++)
                        {
                            Console.WriteLine("k= {0}, x= {1} :", k,x);
                            Console.WriteLine(masr[k]);
                            Console.WriteLine(masg[k]);
                            Console.WriteLine(masb[k]);
                            Console.WriteLine();
                        }
                    for (int j = 0; j < n; j++)
                        for (int i = 0; i < n; i++)
                        {
                            //p.r = BiCubeInterp(masr, x * n + n, y * n + 2 * n, n, x * n + i , y * n + j);
                            p.r = BiCubeInterp(masr, 0, 0, n, i, j);
                            p.b = BiCubeInterp(masb, 0, 0, n, i, j);
                            p.g = BiCubeInterp(masg, 0, 0, n, i, j);
                            /*p.r = BiCubeInterp(masr, 0, 0, n, -i + 2 * n  , j - n);
                            p.b = BiCubeInterp(masb, 0, 0, n, -i + 2 * n , j - n);
                            p.g = BiCubeInterp(masg, 0, 0, n, -i + 2 * n , j - n);*/
                            if (y == 0 && x == 2)
                            {
                                Console.WriteLine("x = {0}  i = {1}  j = {2}", x,i,j);
                                Console.WriteLine(p.r);
                                Console.WriteLine(p.b);
                                Console.WriteLine(p.g);
                            }
                            //p.b = BiCubeInterp(masb, x * n, y * n , n, x * n + i, y * n + j);
                            //p.g = BiCubeInterp(masg, x * n, y * n , n, x * n + i, y * n + j);
                            image1[x * n + j, y * n + i] = p;
                        }
                }
            return image1;

        }

        //this is right version
        static ColorFloatImage BicubeIntInt2(ColorFloatImage image, int n)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width * n, image.Height * n),
                            tmp = new ColorFloatImage(image.Width + 3, image.Height + 3);
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[1 + j, 1 + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, 1 + i] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 2; j++)
                {
                    tmp[j + 1 + image.Width, i + 1] = image[image.Width - 1, i];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 1, i] = image[j, 0];
                }
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 1, i + 1 + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, i] = image[0, 0];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < 2; j++)
                {
                    tmp[j + 1 + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    tmp[j + 1 + image.Width, i + 1 + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, i + 1 + image.Height] = image[0, image.Height - 1];
                }
            float[] masr = new float[16];
            float[] masb = new float[16];
            float[] masg = new float[16];
            for (int y = 0; y < image.Height; y = y + 1)
                for (int x = 0; x < image.Width; x = x + 1)
                {
                    ColorFloatPixel p = image[0, 0];
                    ColorFloatPixel p00 = tmp[x + 1, y + 1];
                    masr[0] = p00.r;
                    masb[0] = p00.b;
                    masg[0] = p00.g;
                    ColorFloatPixel p01 = tmp[x + 1, y + 2];
                    masr[1] = p01.r;
                    masb[1] = p01.b;
                    masg[1] = p01.g;
                    ColorFloatPixel p10 = tmp[x + 2, y + 1];
                    masr[2] = p10.r;
                    masb[2] = p10.b;
                    masg[2] = p10.g;
                    ColorFloatPixel p11 = tmp[x + 2, y + 2];
                    masr[3] = p11.r;
                    masb[3] = p11.b;
                    masg[3] = p11.g;
                    ColorFloatPixel p0_1 = tmp[x + 1, y];
                    masr[4] = p0_1.r;
                    masb[4] = p0_1.b;
                    masg[4] = p0_1.g;
                    ColorFloatPixel p_10 = tmp[x, y + 1];
                    masr[5] = p_10.r;
                    masb[5] = p_10.b;
                    masg[5] = p_10.g;
                    ColorFloatPixel p1_1 = tmp[x + 2, y];
                    masr[6] = p1_1.r;
                    masb[6] = p1_1.b;
                    masg[6] = p1_1.g;
                    ColorFloatPixel p_11 = tmp[x, y + 2];
                    masr[7] = p_11.r;
                    masb[7] = p_11.b;
                    masg[7] = p_11.g;
                    ColorFloatPixel p02 = tmp[x + 1, y + 3];
                    masr[8] = p02.r;
                    masb[8] = p02.b;
                    masg[8] = p02.g;
                    ColorFloatPixel p20 = tmp[x + 3, y + 1];
                    masr[9] = p20.r;
                    masb[9] = p20.b;
                    masg[9] = p20.g;
                    ColorFloatPixel p_1_1 = tmp[x, y];
                    masr[10] = p_1_1.r;
                    masb[10] = p_1_1.b;
                    masg[10] = p_1_1.g;
                    ColorFloatPixel p12 = tmp[x + 2, y + 3];
                    masr[11] = p12.r;
                    masb[11] = p12.b;
                    masg[11] = p12.g;
                    ColorFloatPixel p21 = tmp[x + 3, y + 2];
                    masr[12] = p21.r;
                    masb[12] = p21.b;
                    masg[12] = p21.g;
                    ColorFloatPixel p_12 = tmp[x, y + 3];
                    masr[13] = p_12.r;
                    masb[13] = p_12.b;
                    masg[13] = p_12.g;
                    ColorFloatPixel p2_1 = tmp[x + 3, y];
                    masr[14] = p2_1.r;
                    masb[14] = p2_1.b;
                    masg[14] = p2_1.g;
                    ColorFloatPixel p22 = tmp[x + 3, y + 3];
                    masr[15] = p22.r;
                    masb[15] = p22.b;
                    masg[15] = p22.g;
                    if (y == 0 && x == 2)
                        for (int k = 0; k < 16; k++)
                        {
                            Console.WriteLine("k= {0}, x= {1} :", k, x);
                            Console.WriteLine(masr[k]);
                            Console.WriteLine(masg[k]);
                            Console.WriteLine(masb[k]);
                            Console.WriteLine();
                        }
                    for (int j = 0; j < n; j++)
                        for (int i = 0; i < n; i++)
                        {
                            //p.r = BiCubeInterp(masr, x * n + n, y * n + 2 * n, n, x * n + i , y * n + j);
                            p.r = BiCubeInterp(masr, 0, 0, n, i, j);
                            p.b = BiCubeInterp(masb, 0, 0, n, i, j);
                            p.g = BiCubeInterp(masg, 0, 0, n, i, j);
                            /*p.r = BiCubeInterp(masr, 0, 0, n, -i + 2 * n  , j - n);
                            p.b = BiCubeInterp(masb, 0, 0, n, -i + 2 * n , j - n);
                            p.g = BiCubeInterp(masg, 0, 0, n, -i + 2 * n , j - n);*/
                            if (y == 0 && x == 2)
                            {
                                Console.WriteLine("x = {0}  i = {1}  j = {2}", x, i, j);
                                Console.WriteLine(p.r);
                                Console.WriteLine(p.b);
                                Console.WriteLine(p.g);
                            }
                            //p.b = BiCubeInterp(masb, x * n, y * n , n, x * n + i, y * n + j);
                            //p.g = BiCubeInterp(masg, x * n, y * n , n, x * n + i, y * n + j);
                            image1[x * n + j, y * n + i] = p;
                        }
                }
            return image1;

        }

        static ColorFloatImage BicubeIntInt(ColorFloatImage image, int n)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width * n, image.Height * n),
                            tmp = new ColorFloatImage(image.Width + 5 * n, image.Height + 5 * n);
            /*for (int y = 0; y < image.Height; y++)
                for (int x = 0; x < image.Width; x++)
                {
                    image1[x * n, y * n] = image[x, y];
                }*/
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j * n, i * n] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 8 * n; j++)
                {
                    tmp[j + n*image.Width, i*n] = image[image.Width - 1, i];
                }
            for (int i = 0; i < 8 * n; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j * n, i + n * image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < 8 * n; i++)
                for (int j = 0; j < 8 * n; j++)
                {
                    tmp[j + n * image.Width, i + n * image.Height] = image[image.Width - 1, image.Height - 1];
                }
            //tmp[image.Width * n, image.Height * n] = image[image.Width - 1, image.Height - 1];
            float[] masr = new float[16];
            float[] masb = new float[16];
            float[] masg = new float[16];
            for (int y = 0; y < image.Height+4; y= y+3)
                for (int x = 0; x < image.Width+4; x= x+3)
                {
                    ColorFloatPixel p = image[0, 0];
                    ColorFloatPixel p00 = tmp[x * n + n, y * n + n];
                    masr[0] = p00.r;
                    masb[0] = p00.b;
                    masg[0] = p00.g;
                    ColorFloatPixel p01 = tmp[x * n + n, y * n + 2*n];
                    masr[1] = p01.r;
                    masb[1] = p01.b;
                    masg[1] = p01.g;
                    ColorFloatPixel p10 = tmp[x * n + 2*n, y * n + n];
                    masr[2] = p10.r;
                    masb[2] = p10.b;
                    masg[2] = p10.g;
                    ColorFloatPixel p11 = tmp[x * n + 2 * n, y * n + 2*n];
                    masr[3] = p11.r;
                    masb[3] = p11.b;
                    masg[3] = p11.g;
                    ColorFloatPixel p0_1 = tmp[x * n + n, y * n];
                    masr[4] = p0_1.r;
                    masb[4] = p0_1.b;
                    masg[4] = p0_1.g;
                    ColorFloatPixel p_10 = tmp[x * n, y * n + n];
                    masr[5] = p_10.r;
                    masb[5] = p_10.b;
                    masg[5] = p_10.g;
                    ColorFloatPixel p1_1 = tmp[x * n + 2*n, y * n];
                    masr[6] = p1_1.r;
                    masb[6] = p1_1.b;
                    masg[6] = p1_1.g;
                    ColorFloatPixel p_11 = tmp[x * n, y * n + 2*n];
                    masr[7] = p_11.r;
                    masb[7] = p_11.b;
                    masg[7] = p_11.g;
                    ColorFloatPixel p02 = tmp[x * n + n, y * n + 3*n];
                    masr[8] = p02.r;
                    masb[8] = p02.b;
                    masg[8] = p02.g;
                    ColorFloatPixel p20 = tmp[x * n + 3 * n, y * n];
                    masr[9] = p20.r;
                    masb[9] = p20.b;
                    masg[9] = p20.g;
                    ColorFloatPixel p_1_1 = tmp[x * n , y * n];
                    masr[10] = p_1_1.r;
                    masb[10] = p_1_1.b;
                    masg[10] = p_1_1.g;
                    ColorFloatPixel p12 = tmp[x * n + 2*n, y * n + 3*n];
                    masr[11] = p12.r;
                    masb[11] = p12.b;
                    masg[11] = p12.g;
                    ColorFloatPixel p21 = tmp[x * n + 3*n, y * n + 2*n];
                    masr[12] = p21.r;
                    masb[12] = p21.b;
                    masg[12] = p21.g;
                    ColorFloatPixel p_12 = tmp[x * n, y * n + 3*n];
                    masr[13] = p_12.r;
                    masb[13] = p_12.b;
                    masg[13] = p_12.g;
                    ColorFloatPixel p2_1 = tmp[x * n + 3*n, y * n];
                    masr[14] = p2_1.r;
                    masb[14] = p2_1.b;
                    masg[14] = p2_1.g;
                    ColorFloatPixel p22 = tmp[x * n + 3*n, y * n +3*n];
                    masr[15] = p22.r;
                    masb[15] = p22.b;
                    masg[15] = p22.g;
                    if(y==0 && x==0)
                        for (int k = 0; k < 16; k++)
                        {
                            Console.WriteLine("k= {0} :", k);
                            Console.WriteLine(masr[k]);
                            Console.WriteLine(masg[k]);
                            Console.WriteLine(masb[k]);
                            Console.WriteLine();
                        }
                    for (int j = 0; j < 3*n; j++)
                        for (int i = 0; i < 3*n; i++)
                        {
                            //p.r = BiCubeInterp(masr, x * n + n, y * n + 2 * n, n, x * n + i , y * n + j);
                            p.r = BiCubeInterp(masr, 0, 0, n, i - n, j - n);
                            p.b = BiCubeInterp(masb, 0, 0, n, i - n, j - n);
                            p.g = BiCubeInterp(masg, 0, 0, n, i - n, j - n);
                            /*p.r = BiCubeInterp(masr, 0, 0, n, -i + 2 * n  , j - n);
                            p.b = BiCubeInterp(masb, 0, 0, n, -i + 2 * n , j - n);
                            p.g = BiCubeInterp(masg, 0, 0, n, -i + 2 * n , j - n);*/
                            /*if (y == 0)
                            {
                                Console.WriteLine(p.r);
                                Console.WriteLine(p.b);
                                Console.WriteLine(p.g);
                            }*/
                            //p.b = BiCubeInterp(masb, x * n, y * n , n, x * n + i, y * n + j);
                            //p.g = BiCubeInterp(masg, x * n, y * n , n, x * n + i, y * n + j);
                            tmp[x * n + j, y * n + i] = p;
                            if (x == 0 && y == 0)
                            {
                                Console.WriteLine("after {0} {1}:", i,j);
                                Console.WriteLine(p.r);
                                Console.WriteLine(p.g);
                                Console.WriteLine(p.b);
                                Console.WriteLine();
                            }
                        }
                    /*for (int i = 1; i < n; i++)
                    {
                        p.r = BiLinearInterp(p12.r, p22.r, p11.r, p21.r, x * n, y * n + n, n, x * n + i, y * n);
                        p.b = BiLinearInterp(p12.b, p22.b, p11.b, p21.b, x * n, y * n + n, n, x * n + i, y * n);
                        p.g = BiLinearInterp(p12.g, p22.g, p11.g, p21.g, x * n, y * n + n, n, x * n + i, y * n);
                        tmp[x * n + i, y * n] = p;
                    }
                    for (int i = 1; i < n; i++)
                    {
                        p.r = BiLinearInterp(p12.r, p22.r, p11.r, p21.r, x * n, y * n + n, n, x * n, y * n + i);
                        p.b = BiLinearInterp(p12.b, p22.b, p11.b, p21.b, x * n, y * n + n, n, x * n, y * n + i);
                        p.g = BiLinearInterp(p12.g, p22.g, p11.g, p21.g, x * n, y * n + n, n, x * n, y * n + i);
                        tmp[x * n, y * n + i] = p;
                    }*/
                }
            for (int y = 0; y < image.Height * n; y++)
                for (int x = 0; x < image.Width * n; x++)
                {
                    image1[x, y] = tmp[x, y];
                }
            return image1;
        }

        static ColorFloatImage BicubeIntReal1(ColorFloatImage image, double p)
        {
            ColorFloatImage image1 = new ColorFloatImage((int)(image.Width * p), (int)(image.Height * p)),
                            tmp = new ColorFloatImage(image.Width + 5, image.Height + 5);
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[1 + j, 1 + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, 1 + i] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 4; j++)
                {
                    tmp[j + 1 + image.Width, i + 1] = image[image.Width - 1, i];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 1, i] = image[j, 0];
                }
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 1, i + 1 + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, i] = image[0, 0];
                }
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < 4; j++)
                {
                    tmp[j + 1 + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                {
                    tmp[j + 1 + image.Width, i + 1 + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 1; j++)
                {
                    tmp[j, i + 1 + image.Height] = image[0, image.Height - 1];
                }
            float[] masr = new float[16];
            float[] masb = new float[16];
            float[] masg = new float[16];
            int l = 0, m = 0;
            for (int y = 0; y < (int)(image.Height*p); y = y + 1)
                for (int x = 0; x < (int)(image.Width*p); x = x + 1)
                {
                    l = (int)Math.Floor(y / p);
                    m = (int)Math.Floor(x / p);
                    ColorFloatPixel p1 = image[0, 0];
                    ColorFloatPixel p00 = tmp[m + 1, l + 1];
                    masr[0] = p00.r;
                    masb[0] = p00.b;
                    masg[0] = p00.g;
                    ColorFloatPixel p01 = tmp[m + 1, l + 2];
                    masr[1] = p01.r;
                    masb[1] = p01.b;
                    masg[1] = p01.g;
                    ColorFloatPixel p10 = tmp[m + 2, l + 1];
                    masr[2] = p10.r;
                    masb[2] = p10.b;
                    masg[2] = p10.g;
                    ColorFloatPixel p11 = tmp[m + 2, l + 2];
                    masr[3] = p11.r;
                    masb[3] = p11.b;
                    masg[3] = p11.g;
                    ColorFloatPixel p0_1 = tmp[m + 1, l];
                    masr[4] = p0_1.r;
                    masb[4] = p0_1.b;
                    masg[4] = p0_1.g;
                    ColorFloatPixel p_10 = tmp[m, l + 1];
                    masr[5] = p_10.r;
                    masb[5] = p_10.b;
                    masg[5] = p_10.g;
                    ColorFloatPixel p1_1 = tmp[m + 2, l];
                    masr[6] = p1_1.r;
                    masb[6] = p1_1.b;
                    masg[6] = p1_1.g;
                    ColorFloatPixel p_11 = tmp[m, l + 2];
                    masr[7] = p_11.r;
                    masb[7] = p_11.b;
                    masg[7] = p_11.g;
                    ColorFloatPixel p02 = tmp[m + 1, l + 3];
                    masr[8] = p02.r;
                    masb[8] = p02.b;
                    masg[8] = p02.g;
                    ColorFloatPixel p20 = tmp[m + 3, l + 1];
                    masr[9] = p20.r;
                    masb[9] = p20.b;
                    masg[9] = p20.g;
                    ColorFloatPixel p_1_1 = tmp[m, l];
                    masr[10] = p_1_1.r;
                    masb[10] = p_1_1.b;
                    masg[10] = p_1_1.g;
                    ColorFloatPixel p12 = tmp[m + 2, l + 3];
                    masr[11] = p12.r;
                    masb[11] = p12.b;
                    masg[11] = p12.g;
                    ColorFloatPixel p21 = tmp[m + 3, l + 2];
                    masr[12] = p21.r;
                    masb[12] = p21.b;
                    masg[12] = p21.g;
                    ColorFloatPixel p_12 = tmp[m, l + 3];
                    masr[13] = p_12.r;
                    masb[13] = p_12.b;
                    masg[13] = p_12.g;
                    ColorFloatPixel p2_1 = tmp[m + 3, l];
                    masr[14] = p2_1.r;
                    masb[14] = p2_1.b;
                    masg[14] = p2_1.g;
                    ColorFloatPixel p22 = tmp[m + 3, l + 3];
                    masr[15] = p22.r;
                    masb[15] = p22.b;
                    masg[15] = p22.g;
                    if (y == 0 && x == 120)
                        for (int k = 0; k < 16; k++)
                        {
                            Console.WriteLine("k= {0}, x= {1} :", k, x);
                            Console.WriteLine(masr[k]);
                            Console.WriteLine(masg[k]);
                            Console.WriteLine(masb[k]);
                            Console.WriteLine();
                        }
                    if (y == 0 && x == 120)
                        Console.WriteLine("(float)(x/p)={0} , (float)(y/p)={1}", (float)(x / p), (float)(y / p));
                            //p.r = BiCubeInterp(masr, x * n + n, y * n + 2 * n, n, x * n + i , y * n + j);
                    p1.r = BiCubeInterp(masr, 0, 0, 1, (float)(y / p) - (float)Math.Floor(y / p) + 1.0f, (float)(x / p) - (float)Math.Floor(x / p) + 1.0f);
                    p1.b = BiCubeInterp(masb, 0, 0, 1, (float)(y / p) - (float)Math.Floor(y / p) + 1.0f, (float)(x / p) - (float)Math.Floor(x / p) + 1.0f);
                    p1.g = BiCubeInterp(masg, 0, 0, 1, (float)(y / p) - (float)Math.Floor(y / p) + 1.0f, (float)(x / p) - (float)Math.Floor(x / p) + 1.0f);
                            /*p.r = BiCubeInterp(masr, 0, 0, n, -i + 2 * n  , j - n);
                            p.b = BiCubeInterp(masb, 0, 0, n, -i + 2 * n , j - n);
                            p.g = BiCubeInterp(masg, 0, 0, n, -i + 2 * n , j - n);*/
                            if (y == 0 && x == 120)
                            {
                                Console.WriteLine("x = {0}  i = {1}  j = {2} floor= {3} ceiling={4}", x, x/p, y/p, (int)Math.Floor(x/p), (int)Math.Ceiling(x/p));
                                Console.WriteLine(p1.r);
                                Console.WriteLine(p1.b);
                                Console.WriteLine(p1.g);
                            }
                            //p.b = BiCubeInterp(masb, x * n, y * n , n, x * n + i, y * n + j);
                            //p.g = BiCubeInterp(masg, x * n, y * n , n, x * n + i, y * n + j);
                    image1[x,y] = p1;
                }
            return image1;

        }
       
        static ColorFloatImage Lanczoh1(ColorFloatImage image, int n)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width * n, image.Height * n), tmp = new ColorFloatImage(image.Width + 2 * 3, image.Height + 2 * 3);
            float[ ] masLx = new float[7];
            float[ ] masLy = new float[7];
            float red=0,green=0, blue=0;
            ColorFloatPixel p=image[0,0];
            for (int i = -3; i<4; i++)
            {
                if(i==0)
                {
                    masLx[i+3]=1;
                    masLy[i+3]=1;
                }
                else
                {
                    masLx[i+3]=(float)((3*Math.Sin(Math.PI*i)*Math.Sin(Math.PI*i/3))/((Math.PI)*(Math.PI)*(i)*(i)));
                    masLy[i+3]=(float)((3*Math.Sin(Math.PI*i)*Math.Sin(Math.PI*i/3))/((Math.PI)*(Math.PI)*(i)*(i)));
                }
            }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[3 + j, 3 + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j, i + 3] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j + 3 + image.Width, i + 3] = image[image.Width - 1, i];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 3, i] = image[j, 0];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 3, i + 3 + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j, i] = image[0, 0];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j + 3 + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j + 3 + image.Width, i + 3 + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j, i + 3 + image.Height] = image[0, image.Height - 1];
                }
            for (int y=0; y<image.Height*n; y++)
                for (int x = 0; x < image.Width*n; x++)
                {
                    if (x % n == 0)
                    {
                        if (y % n == 0)
                            image1[x, y] = tmp[x / n + 3, y / n + 3];
                        else
                        {
                            red=blue=green=0;
                            for(int i=-3;i<4;i++)
                                for(int j=-3;j<4;j++)
                                {
                                    red = red + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].r * masLy[i + 3];
                                    blue = blue + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].b * masLy[i + 3];
                                    green = green + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].g * masLy[i + 3];
                                }
                            p.r=red;
                            p.b=blue;
                            p.g=green;
                            image1[x, y]=p;
                        }
                    }
                    else
                    {
                        red = blue = green = 0;
                        for (int i = -3; i < 4; i++)
                            for (int j = -3; j < 4; j++)
                            {
                                red = red + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].r * masLy[i + 3];
                                blue = blue + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].b * masLy[i + 3];
                                green = green + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].g * masLy[i + 3];
                            }
                        p.r = red;
                        p.b = blue;
                        p.g = green;
                        image1[x, y] = p;
                    }
                }
            return image1;
        }

        static ColorFloatImage Lanczoh(ColorFloatImage image, double q)
        {
            ColorFloatImage image1 = new ColorFloatImage((int)(image.Width * q), (int)(image.Height * q)), 
                            tmp = new ColorFloatImage(image.Width + 2 * 3, image.Height + 2 * 3);
            float[] masLx = new float[7];
            float[] masLy = new float[7];
            float red = 0, green = 0, blue = 0;
            ColorFloatPixel p = image[0, 0];
            /*for (int i = -3; i < 4; i++)
            {
                if (i == 0)
                {
                    masLx[i + 3] = 1;
                    masLy[i + 3] = 1;
                }
                else
                {
                    masLx[i + 3] = (float)((3 * Math.Sin(Math.PI * i) * Math.Sin(Math.PI * i / 3)) / ((Math.PI) * (Math.PI) * (i) * (i)));
                    masLy[i + 3] = (float)((3 * Math.Sin(Math.PI * i) * Math.Sin(Math.PI * i / 3)) / ((Math.PI) * (Math.PI) * (i) * (i)));
                }
            }*/
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[3 + j, 3 + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j, i + 3] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j + 3 + image.Width, i + 3] = image[image.Width - 1, i];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 3, i] = image[j, 0];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j + 3, i + 3 + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j, i] = image[0, 0];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j + 3 + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j + 3 + image.Width, i + 3 + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    tmp[j, i + 3 + image.Height] = image[0, image.Height - 1];
                }
            for (int y = 0; y < (int)(image.Height * q); y++)
                for (int x = 0; x < (int)(image.Width * q); x++)
                {
                    /*if (x % n == 0)
                    {
                        if (y % n == 0)
                            image1[x, y] = tmp[x / n + 3, y / n + 3];
                        else
                        {*/
                    for (int i = -3; i < 4; i++)
                    {
                        if ((x / q - Math.Floor(x / q) + i) == 0)
                        {
                            masLx[i + 3] = 1;
                            //masLy[i + 3] = 1;
                        }
                        else
                        {
                            masLx[i + 3] = (float)((3 * Math.Sin(Math.PI * (x / q - Math.Floor(x / q) + i)) 
                                           * Math.Sin(Math.PI * (x / q - Math.Floor(x / q) + i) / 3))
                                           / ((Math.PI) * (Math.PI) * (x / q - Math.Floor(x / q) + i) * (x / q - Math.Floor(x / q) + i)));
                            //masLy[i + 3] = (float)((3 * Math.Sin(Math.PI * (y / q - Math.Floor(y / q) + i)) * Math.Sin(Math.PI * (y / q - Math.Floor(y / q) + i) / 3))
                                           /// ((Math.PI) * (Math.PI) * (y / q - Math.Floor(y / q) + i) * (y / q - Math.Floor(y / q) + i)));
                        }
                        if ((y / q - Math.Floor(y / q) + i) == 0)
                        {
                            masLy[i + 3] = 1;
                        }
                        else
                        {
                            masLy[i + 3] = (float)((3 * Math.Sin(Math.PI * (y / q - Math.Floor(y / q) + i))
                                           * Math.Sin(Math.PI * (y / q - Math.Floor(y / q) + i) / 3))
                                           / ((Math.PI) * (Math.PI) * (y / q - Math.Floor(y / q) + i) * (y / q - Math.Floor(y / q) + i)));
                        }
                    } 
                            if(x==2 && y==4)
                                for (int i = 0; i < 7; i++)
                                {
                                    Console.WriteLine("masly= {0}, maslx= {1}",masLx[i],masLy[i]);
                                }
                            red = blue = green = 0;
                            for (int i = -3; i < 4; i++)
                                for (int j = -3; j < 4; j++)
                                {
                                    red = red + masLx[j + 3] * tmp[(int)(Math.Floor(x / q)) - j + 3, (int)(Math.Floor(y / q)) - i + 3].r * masLy[i + 3];
                                    blue = blue + masLx[j + 3] * tmp[(int)(Math.Floor(x / q)) - j + 3, (int)(Math.Floor(y / q)) - i + 3].b * masLy[i + 3];
                                    green = green + masLx[j + 3] * tmp[(int)(Math.Floor(x / q)) - j + 3, (int)(Math.Floor(y / q)) - i + 3].g * masLy[i + 3];
                                }
                            p.r = red;
                            p.b = blue;
                            p.g = green;
                            image1[x, y] = p;
                      /*  }
                    }
                    else
                    {
                        red = blue = green = 0;
                        for (int i = -3; i < 4; i++)
                            for (int j = -3; j < 4; j++)
                            {
                                red = red + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].r * masLy[i + 3];
                                blue = blue + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].b * masLy[i + 3];
                                green = green + masLx[j + 3] * tmp[x / n + j + 3, y / n + i + 3].g * masLy[i + 3];
                            }
                        p.r = red;
                        p.b = blue;
                        p.g = green;
                        image1[x, y] = p;
                    }*/
                }
            return image1;
        }

        static ColorFloatImage LowInt(ColorFloatImage image, int n)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width / n, image.Height / n);
            image= Gauss(image, (int) (Math.Sqrt(n*n-1)+1),Math.Sqrt(n*n-1));
            for (int i = 0; i < image.Height / n; i++)
                for (int j = 0; j < image.Width / n; j++)
                {
                    image1[j, i] = image[j*n, i*n];
                }
            return image1;
        }

        static float BiLinearInterp1(float y12, float y22, float y11, float y21, int x1, int y1, int n, double x, double y)
        {
            return (float)((-y11) * (x1 + n - x) * (y1 - n - y) + (-y21) * (x - x1) * (y1 - n - y) + 
                (-y12) * (x1 + n - x) * (y - y1) + (-y22) * (x - x1) * (y - y1)) / (n * n);
        }
        
        static ColorFloatImage LowDouble(ColorFloatImage image, double p)
        {
            ColorFloatImage image1 = new ColorFloatImage((int) (image.Width / p), (int) (image.Height / p)),
                            tmp = new ColorFloatImage(image.Width + 1, image.Height + 1); 
            image = Gauss(image, (int)(Math.Sqrt(p * p - 1) + 1), Math.Sqrt(p * p - 1));
            ColorFloatPixel p1 = image[0, 0],
                p11 = image[0, 0],
                p12 = image[0, 0],
                p21 = image[0, 0],
                p22 = image[0, 0];
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j, i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
            {
                tmp[image.Width , i] = image[image.Width - 1, i];
            }
            for (int j = 0; j < image.Width; j++)
            {
                tmp[j, image.Height] = image[j, image.Height - 1];
            }
            tmp[image.Width, image.Height ] = image[image.Width - 1, image.Height - 1];
                 
            for (int i = 0; i < (int)(image.Height / p); i++)
                for (int j = 0; j < (int)(image.Width / p); j++)
                {
                    //image1[j, i] = image[j * n, i * n];
                    p11 = tmp[(int)Math.Floor(j * p), (int)(Math.Ceiling(i * p))];
                    p12 = tmp[(int)(Math.Floor(j * p)), (int)(Math.Floor(i * p))];
                    p21 = tmp[(int)(Math.Ceiling(j * p)), (int)(Math.Ceiling(i * p))];
                    p22 = tmp[(int)(Math.Ceiling(j * p)), (int)(Math.Floor(i * p))];
                    p1.r = BiLinearInterp1(p12.r, p22.r, p11.r, p21.r, (int)(Math.Floor(j * p)), (int)(Math.Ceiling(i * p)), 1, j * p, i * p);
                    p1.b = BiLinearInterp1(p12.b, p22.b, p11.b, p21.b, (int)(Math.Floor(j * p)), (int)(Math.Ceiling(i * p)), 1, j * p, i * p);
                    p1.g = BiLinearInterp1(p12.g, p22.g, p11.g, p21.g, (int)(Math.Floor(j * p)), (int)(Math.Ceiling(i * p)), 1, j * p, i * p);
                    image1[j, i] = p1;
                }
            return image1;
        }

        static ColorFloatImage BILinearReal(ColorFloatImage image, double p)
        {
            ColorFloatImage image1 = new ColorFloatImage((int)(image.Width * p), (int)(image.Height * p)),
                            tmp = new ColorFloatImage(image.Width + 1, image.Height + 1);
            ColorFloatPixel p1 = image[0, 0],
                p11 = image[0, 0],
                p12 = image[0, 0],
                p21 = image[0, 0],
                p22 = image[0, 0];
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmp[j, i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
            {
                tmp[image.Width, i] = image[image.Width - 1, i];
            }
            for (int j = 0; j < image.Width; j++)
            {
                tmp[j, image.Height] = image[j, image.Height - 1];
            }
            tmp[image.Width, image.Height] = image[image.Width - 1, image.Height - 1];
            for (int i = 0; i < (int)(image.Height * p); i++)
                for (int j = 0; j < (int)(image.Width * p); j++)
                {
                    //image1[j, i] = image[j * n, i * n];
                    p11 = tmp[(int)Math.Floor(j / p), (int)(Math.Ceiling(i / p))];
                    p12 = tmp[(int)(Math.Floor(j / p)), (int)(Math.Floor(i / p))];
                    p21 = tmp[(int)(Math.Ceiling(j / p)), (int)(Math.Ceiling(i / p))];
                    p22 = tmp[(int)(Math.Ceiling(j / p)), (int)(Math.Floor(i / p))];
                    p1.r = BiLinearInterp1(p12.r, p22.r, p11.r, p21.r, (int)(Math.Floor(j / p)), (int)(Math.Ceiling(i / p)), 1, j / p, i / p);
                    p1.b = BiLinearInterp1(p12.b, p22.b, p11.b, p21.b, (int)(Math.Floor(j / p)), (int)(Math.Ceiling(i / p)), 1, j / p, i / p);
                    p1.g = BiLinearInterp1(p12.g, p22.g, p11.g, p21.g, (int)(Math.Floor(j / p)), (int)(Math.Ceiling(i / p)), 1, j / p, i / p);
                    image1[j, i] = p1;
                }
            return image1;
        }

        static ColorFloatImage GaussGradNMS(ColorFloatImage image, int n, double sigma)
        {
            ColorFloatImage image1 = new ColorFloatImage(image.Width, image.Height), tmpx = new ColorFloatImage(image.Width + 2 * n, image.Height + 2 * n),
            tmpy = new ColorFloatImage(image.Width + 2 * n, image.Height + 2 * n);
            float[] masGx = new float[(2 * n + 1)];
            int k1 = 2 * n + 1;
            float coeff = (float)(-1.0f / (2.0f * Math.PI * sigma * sigma * sigma * sigma));
            float sumx = 0;
            for (int i = 0; i < 2 * n + 1; i++)
            {
                masGx[i] = coeff * (i - n) * (float)Math.Exp(-(Math.Pow(n - i, 2)) / (2 * sigma * sigma));
                sumx += Math.Abs(masGx[i]);
                Console.WriteLine(masGx[i]);
            }
            for (int i = 0; i < 2 * n + 1; i++)
            {
                masGx[i] = masGx[i] / sumx;

            }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmpx[n + j, n + i] = image[j, i];
                    tmpy[n + j, n + i] = image[j, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < n; j++)
                {
                    tmpx[j, n + i] = image[0, i];
                    tmpy[j, n + i] = image[0, i];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < n; j++)
                {
                    tmpx[j + n + image.Width, n + i] = image[image.Width - 1, i];
                    tmpy[j + n + image.Width, n + i] = image[image.Width - 1, i];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmpx[j + n, i] = image[j, 0];
                    tmpy[j + n, i] = image[j, 0];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    tmpx[j + n, i + n + image.Height] = image[j, image.Height - 1];
                    tmpy[j + n, i + n + image.Height] = image[j, image.Height - 1];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmpx[j, i] = image[0, 0];
                    tmpy[j, i] = image[0, 0];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmpx[j + n + image.Width, i] = image[image.Width - 1, 0];
                    tmpy[j + n + image.Width, i] = image[image.Width - 1, 0];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmpx[j + n + image.Width, i + n + image.Height] = image[image.Width - 1, image.Height - 1];
                    tmpy[j + n + image.Width, i + n + image.Height] = image[image.Width - 1, image.Height - 1];
                }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    tmpx[j, i + n + image.Height] = image[0, image.Height - 1];
                    tmpy[j, i + n + image.Height] = image[0, image.Height - 1];
                }
            for (int i = 0; i < image.Height; i++)
                for (int j = 0; j < image.Width; j++)
                {
                    ColorFloatPixel px = image[j, i], py = image[j, i], p = image[j, i];
                    p.r = py.r = px.r = 0;
                    p.g = py.g = px.g = 0;
                    p.b = py.b = px.b = 0;
                    for (int k = 0; k < 2 * n + 1; k++)
                    //for (int m = 0; m < 2 * n + 1; m++)
                    {
                        ColorFloatPixel p1x = tmpx[k + j, n + i], p1y = tmpy[n + j, k + i];
                        px.r = px.r + p1x.r * masGx[k];
                        px.g = px.g + p1x.g * masGx[k];
                        px.b = px.b + p1x.b * masGx[k];
                        py.r = py.r + p1y.r * masGx[k];
                        py.g = py.g + p1y.g * masGx[k];
                        py.b = py.b + p1y.b * masGx[k];
                    }
                    p.r = (float)Math.Sqrt(px.r * px.r + py.r * py.r);
                    p.b = (float)Math.Sqrt(px.b * px.b + py.b * py.b);
                    p.g = (float)Math.Sqrt(px.g * px.g + py.g * py.g);
                    image1[j, i] = p;
                }
            double gx = 0, gy = 0, thetta=0, m0=0 , m1=0, m2=0;
            for (int i = 1; i < image.Height-1; i++)
                for (int j = 1; j < image.Width-1; j++)
                {
                    ColorFloatPixel px = image[j, i], py = image[j, i], p = image[j, i];
                    p.r = py.r = px.r = 0;
                    p.g = py.g = px.g = 0;
                    p.b = py.b = px.b = 0;
                    for (int k = 0; k < 2 * n + 1; k++)
                    //for (int m = 0; m < 2 * n + 1; m++)
                    {
                        ColorFloatPixel p1x = tmpx[k + j, n + i], p1y = tmpy[n + j, k + i];
                        px.r = px.r + p1x.r * masGx[k];
                        px.g = px.g + p1x.g * masGx[k];
                        px.b = px.b + p1x.b * masGx[k];
                        py.r = py.r + p1y.r * masGx[k];
                        py.g = py.g + p1y.g * masGx[k];
                        py.b = py.b + p1y.b * masGx[k];
                    }
                    gx = Math.Sqrt(px.r * px.r + px.b * px.b + px.g * px.g);
                    gy = Math.Sqrt(py.r * py.r + py.b * py.b + py.g * py.g);
                    thetta = Math.Atan2(gy, gx);
                    if ((thetta > 0 && thetta <= (Math.PI) / 8) || ((thetta > 7 * Math.PI / 8) && (thetta <= Math.PI)))
                    {
                        p = image1[j - 1, i];
                        m0 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                        p = image1[j , i];
                        m1 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                        p = image1[j + 1, i];
                        m2 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                        if (m1 < m0 || m1 < m2)
                        {
                            p.r = p.b = p.g = 0;
                            image1[j, i] = p;
                        }
                    }
                    else
                        if (thetta > Math.PI / 8 && thetta <= 3 * Math.PI / 8)
                        {
                            p = image1[j - 1, i + 1];
                            m0 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                            p = image1[j, i];
                            m1 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                            p = image1[j + 1, i -1];
                            m2 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                            if (m1 < m0 || m1 < m2)
                            {
                                p.r = p.b = p.g = 0;
                                image1[j, i] = p;
                            }
                        }
                        else
                            if (thetta > 3 * Math.PI / 8 && thetta <= 5 * Math.PI / 8)
                            {
                                p = image1[j, i - 1];
                                m0 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                                p = image1[j, i];
                                m1 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                                p = image1[j, i + 1];
                                m2 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                                if (m1 < m0 || m1 < m2)
                                {
                                    p.r = p.b = p.g = 0;
                                    image1[j, i] = p;
                                }
                            }
                            else
                                if (thetta > 5 * Math.PI / 8 && thetta <= 7 * Math.PI / 8)
                                {
                                    p = image1[j - 1, i - 1];
                                    m0 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                                    p = image1[j, i];
                                    m1 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                                    p = image1[j + 1, i + 1];
                                    m2 = Math.Sqrt(p.r * p.r + p.b * p.b + p.g * p.g);
                                    if (m1 < m0 || m1 < m2)
                                    {
                                        p.r = p.b = p.g = 0;
                                        image1[j, i] = p;
                                    }
                                }
                }
            return image1;
        }

        static double MSE(ColorFloatImage image1, ColorFloatImage image2)
        {
            ColorFloatPixel p1 = image1[0, 0], p2 = image2[0, 0];
            double res=0;
            for(int i=0; i<image1.Height; i++)
                for (int j = 0; j < image1.Width; j++)
                {
                    p1=image1[j,i];
                    p2=image2[j,i];
                    res += Math.Sqrt((p2.r - p1.r) * (p2.r - p1.r) + (p2.b - p1.b) * (p2.b - p1.b) + (p2.g - p1.g) * (p2.g - p1.g));
                }
            return res / (image1.Height * image1.Width);
        }

        static double PSNR(ColorFloatImage image1, ColorFloatImage image2)
        {
            double q = MSE(image1, image2);
            if (q <= 0)
                return 0;
            else
                return 10 * Math.Log10((Math.Pow(2, 8) - 1 )* (Math.Pow(2, 8)-1) / MSE(image1, image2));
        }

        static double MSSIM(ColorFloatImage image1, ColorFloatImage image2)
        {
            double mux = 0, muy = 0, sigmax = 0, sigmay = 0, sigmaxy = 0;
            double ssim=0,c1 =(0.01*(Math.Pow(2,8)-1))*(0.01*(Math.Pow(2,8)-1)),c2=(0.03*(Math.Pow(2,8)-1))*(0.03*(Math.Pow(2,8)-1));
            ColorFloatPixel p1 = image1[0, 0], p2 = image2[0, 0];
            for(int i=0; i<image1.Height-8; i++)
                for (int j = 0; j < image1.Width - 8; j++)
                {
                    //Console.WriteLine("j={0}; i={1} ", j, i); 
                    mux = 0;
                    muy = 0;
                    sigmax = 0;
                    sigmay = 0;
                    sigmaxy = 0;
                    for(int y=0; y<8; y++)
                        for (int x = 0; x < 8; x++)
                        {
                            p1 = image1[j + x, i + y];
                            mux = mux + Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g)/64;
                            p2 = image2[j + x, i + y];
                            muy = muy + Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g)/64;
                        }
                    //mux = mux / 64;
                    //muy = muy / 64;
                    for (int y = 0; y < 8; y++)
                        for (int x = 0; x < 8; x++)
                        {
                            p1 = image1[j + x, i + y];
                            sigmax = sigmax + (Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) - mux) * (Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) - mux)
                                                                                                                                                                / 64;
                            p2 = image2[j + x, i + y];
                            sigmay = sigmay + (Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) - muy) * (Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) - muy)
                                                                                                                                                                / 64;
                            sigmaxy = sigmaxy + (Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) - mux) * (Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) - muy)
                                                                                                                                                                / 64;
                        }
                    //sigmax = sigmax / 64;
                    //sigmay = sigmay / 64;
                    //sigmaxy = sigmaxy / 64;
                    ssim = ssim + ((2 * mux * muy + c1) * (2 * sigmaxy + c2) / ((mux * mux + muy * muy + c1) * (sigmax + sigmay + c2) 
                                                                                                    * (image1.Height - 8) * (image1.Width - 8)));
                    /*Console.WriteLine("1= {0}; 2={1}; 3= {2}", (2 * mux * muy + c1) * (2 * sigmaxy + c2),
                        ((2 * mux * muy + c1) * (2 * sigmaxy + c2) / ((mux * mux + muy * muy + c1) * (sigmax + sigmay + c2)
                                                                                                    * (image1.Height - 8) * (image1.Width - 8))),
                        ((mux * mux + muy * muy + c1) * (sigmax + sigmay + c2) * (image1.Height - 8) * (image1.Width - 8)));*/
                }
            return ssim;
        }

        static double SSIM(ColorFloatImage image1, ColorFloatImage image2)
        {
            double mux = 0, muy = 0, sigmax = 0, sigmay = 0, sigmaxy = 0;
            double c1 = (0.01 * (Math.Pow(2, 8) - 1)) * (0.01 * (Math.Pow(2, 8) - 1)), c2 = (0.03 * (Math.Pow(2, 8) - 1)) * (0.03 * (Math.Pow(2, 8) - 1));
            ColorFloatPixel p1 = image1[0, 0], p2 = image2[0, 0];
            mux = 0;
            muy = 0;
            sigmax = 0;
            sigmay = 0;
            sigmaxy = 0;
            for (int i = 0; i < image1.Height; i++)
                for (int j = 0; j < image1.Width; j++)
                {
                     p1 = image1[j, i];
                     mux = mux + Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) / (image1.Height * image1.Width);
                     p2 = image2[j, i];
                     muy = muy + Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) / (image1.Height * image1.Width);
                }
            //mux = mux / (image1.Height*image1.Width);
            //muy = muy / (image1.Height*image1.Width);
            for (int i = 0; i < image1.Height; i++)
                for (int j = 0; j < image1.Width; j++)
                {
                     p1 = image1[j, i];
                     sigmax = sigmax + (Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) - mux) * (Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) - mux) 
                                                                                                                              / (image1.Height * image1.Width);
                     p2 = image2[j, i];
                     sigmay = sigmay + (Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) - muy) * (Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) - muy)
                                                                                                                              / (image1.Height * image1.Width);
                     sigmaxy += (Math.Sqrt(p1.r * p1.r + p1.b * p1.b + p1.g * p1.g) - mux) * (Math.Sqrt(p2.r * p2.r + p2.b * p2.b + p2.g * p2.g) - muy)
                                                                                                                              / (image1.Height * image1.Width);
                }
            //sigmax = sigmax / (image1.Height * image1.Width);
            //sigmay = sigmay / (image1.Height * image1.Width);
            //sigmaxy = sigmaxy / (image1.Height * image1.Width);
            return (2 * mux * muy + c1) * (2 * Math.Sqrt(sigmaxy) + c2) / ((mux * mux + muy * muy + c1) * (sigmax + sigmay + c2));
        }

        static void Main(string[] args)
        {
            //ColorFloatImage image1 = ImageIO.FileToColorFloatImage("C:/Users/HP/Downloads/img/bikes2_87_1.bmp");
            //ColorFloatImage image2 = ImageIO.FileToColorFloatImage("C:/Users/HP/Downloads/img/bikes8.bmp");
            Console.WriteLine("Insert full path to the image file (.bmp)");
            //changed!!!!!!!!!!!!!!!!!
            string InputFileName = "C:/Users/HP/Downloads/img_bicubic/bikes_bicubic_before.bmp";//Console.ReadLine();
            if (!File.Exists(InputFileName))
                return;
            Console.WriteLine("Insert full path to the RESULT (or the SECOND image, if you want to compare) image file (.bmp)");
            //changed!!!!!!!!!!!!!!!!!
            string OutputFileName = "C:/Users/HP/Downloads/img_bicubic/bikes_bicubic.bmp";//Console.ReadLine();
            Console.WriteLine("Insert exact word of the method, which you want to conduct with the image:");
            Console.WriteLine(" BiLinearInt BiCubeInt BilLinearReal BiCubeReal LowInt LowReal Lancozh MSE PSNR SSIM MSSIM Non_Maximum_Supression");
            //changed!!!!!!!!!!!!!!!!1
            string Method = "PSNR";//Console.ReadLine();

            //ColorFloatImage image = ImageIO.FileToColorFloatImage("C:/Users/HP/Downloads/img/Mae4r7Pmnn4.bmp");
            ColorFloatImage image = ImageIO.FileToColorFloatImage(InputFileName);
            if (Method == "BiLinearInt")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" integer value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    int m1 = Convert.ToInt32(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = BileanIntInt(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "BiCubeInt")
            {
                //try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" integer value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    int m1 = Convert.ToInt32(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = BicubeIntInt2(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                //catch (Exception e)
                //{
                //    Console.WriteLine(e.Message);
                //    return;
                //}
            }
            if (Method == "BiLinearReal")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" double value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    double m1 = Convert.ToDouble(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = BILinearReal(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "BiCubeReal")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" double value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    double m1 = Convert.ToDouble(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = BicubeIntReal1(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "Non_Maximum_Supression")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" Int value (sigma>0), the decimal separator is: ','!!!");
                    string Params = Console.ReadLine();
                    double m1 = Convert.ToDouble(Params);
                    if (m1 <= 0)
                        return;
                    ColorFloatImage image1 = GaussGradNMS(image, 3, 1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    return;
                }
            }
            if (Method == "LowInt")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" integer value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    int m1 = Convert.ToInt32(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = LowInt(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "LowReal")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" double value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    double m1 = Convert.ToDouble(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = LowDouble(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "Lancozh")
            {
                try
                {
                    Console.WriteLine("Insert exact word of the parameter, which you want to use:");
                    Console.WriteLine(" double value (>1)!!!");
                    string Params = Console.ReadLine();
                    Console.WriteLine(Params);
                    double m1 = Convert.ToDouble(Params);
                    if (m1 < 1)
                        return;
                    ColorFloatImage image1 = Lanczoh(image,m1);
                    ImageIO.ImageToFile(image1, OutputFileName);
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "MSE")
            {
                try
                {
                    if (!File.Exists(OutputFileName))
                        return;
                     ColorFloatImage image1 = ImageIO.FileToColorFloatImage(OutputFileName);
                    if(image.Width != image1.Width || image.Height !=image1.Height)
                        return;
                    Console.WriteLine("Result = ");
                    Console.WriteLine(MSE(image, image1));
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "PSNR")
            {
                try
                {
                    if (!File.Exists(OutputFileName))
                        return;
                    ColorFloatImage image1 = ImageIO.FileToColorFloatImage(OutputFileName);
                    if (image.Width != image1.Width || image.Height != image1.Height)
                        return;
                    Console.WriteLine("Result = ");
                    Console.WriteLine(PSNR(image, image1));
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "SSIM")
            {
                try
                {
                    if (!File.Exists(OutputFileName))
                        return;
                    ColorFloatImage image1 = ImageIO.FileToColorFloatImage(OutputFileName);
                    if (image.Width != image1.Width || image.Height != image1.Height)
                        return;
                    Console.WriteLine("Result = ");
                    Console.WriteLine(SSIM(image, image1));
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            if (Method == "MSSIM")
            {
                try
                {
                    if (!File.Exists(OutputFileName))
                        return;
                    ColorFloatImage image1 = ImageIO.FileToColorFloatImage(OutputFileName);
                    if (image.Width != image1.Width || image.Height != image1.Height)
                        return;
                    Console.WriteLine("Result = ");
                    Console.WriteLine(MSSIM(image, image1));
                }
                catch (Exception)
                {
                    Console.WriteLine("exception");
                    return;
                }
            }
            //Rotate180Image(image);
            //ColorFloatImage image3 = LowDouble(image1,2.87);
            //ColorFloatImage image2 = SobelVert(image);
            //int sigma=2;
            //ColorFloatImage image3 = GaussGrad(image,sigma,sigma);
            //ImageIO.ImageToFile(image3, "C:/Users/HP/Downloads/img/bike.bmp");
            //Console.WriteLine(MSSIM(image1,image2));
            //Console.ReadLine();
            //ImageIO.ImageToFile(image1, "C:/Users/HP/Downloads/img/bikes1.bmp");
            //ColorFloatImage image = ImageIO.FileToColorFloatImage("C:/Users/HP/Downloads/img/bikes.bmp");

        }
    }
}
