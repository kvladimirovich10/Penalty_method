package main;

import java.util.ArrayList;

import static java.lang.Math.max;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.lang.StrictMath.*;

public class Main {

    private static boolean result = false;

    public static double fun(Matrix x) {


        return pow(x.getByIndex(0, 0), 3) +
                2 * pow(x.getByIndex(1, 0), 2) * x.getByIndex(2, 0) +
                2 * x.getByIndex(2, 0);


    }

    public static double alpha(Matrix x) {


        return pow(max(0, pow(x.getByIndex(0, 0), 2) - x.getByIndex(1, 0) + 2 * x.getByIndex(2, 0) - 2), 2) +
                pow(max(0, pow(x.getByIndex(0, 0), 2) + pow(x.getByIndex(1, 0), 2) + pow(x.getByIndex(2, 0), 2) - 5), 2) +
                pow(abs(pow(x.getByIndex(0, 0), 2) + x.getByIndex(1, 0) + pow(x.getByIndex(2, 0), 2)), 2);

    }

    public static double Falpha(Matrix x, double m) {
        return fun(x) + m * alpha(x);
    }

    public static Matrix CoordSteps(double[] lambda, Matrix[] s, double alpha, double beta, int n, Matrix y, double[] sumStep, boolean[] history, Matrix yResult, double m) {
        int counter = 0;

        yResult.copy(y);

        for (int j = 0 ; j < n ; j++) {
            boolean oldHistory = history[j];

            Matrix op = yResult.plus(s[j].mult(lambda[j]));

            if (Falpha(op, m) < Falpha(yResult, m)) {
                sumStep[j] += lambda[j];
                yResult = op;
                lambda[j] *= alpha;
                history[j] = true;
            } else {
                lambda[j] *= beta;
                history[j] = false;
            }

            if (oldHistory == true && history[j] == false) {
                counter++;
            } else
                counter--;
        }

        if (counter == n) {
            //после каждой удачи последовала неудача
            result = true;
        }
        result = false;

        return yResult;
    }

    public static Matrix MR1(int n, Matrix x0, double epsilon, double m) {
        Matrix[] s = new Matrix[n];

        for (int i = 0 ; i < n ; i++) {
            s[i] = Matrix.getIdentityVector(n, i);
        }

        double[] lambda = new double[n];

        for (int i = 0 ; i < n ; i++) {
            lambda[i] = 0.1;
        }

        Matrix xk = null, xk1 = null;
        double sumStep[] = new double[n];
        boolean history[] = new boolean[n];
        double alpha = 3, beta = -0.5;
        Matrix y1 = x0;
        Matrix yResult = y1.getCopy();

        xk = x0;

        while (true) {
            boolean result = false;
            while (!result) {
                yResult = CoordSteps(lambda, s, alpha, beta, n, y1, sumStep, history, yResult, m);
                if (Falpha(yResult, m) < Falpha(y1, m)) {
                    y1 = yResult;
                    continue;
                }
                if (Falpha(yResult, m) == Falpha(y1, m)) {
                    if (Falpha(yResult, m) < Falpha(xk, m)) {
                        xk1 = yResult;
                        if (norm(xk.minus(xk1), n) < epsilon) {
                            //   System.out.println("x = [ " + xk1.getByIndex(0, 0) + " , " + xk1.getByIndex(1, 0) + "]");
                            return xk1;
                        } else {
                            xk = xk1;
                            break;
                        }
                    }
                }
                if (Falpha(yResult, m) == Falpha(xk, m)) {
                    boolean isExit = true;
                    for (double current : lambda) {
                        if (abs(current) > epsilon) {
                            isExit = false;
                        }
                    }
                    if (isExit) {
                        // System.out.println("x = [ " + xk.getByIndex(0, 0) + " , " + xk.getByIndex(1, 0) + "]");
                        return xk;
                    } else {
                        y1 = yResult;
                        continue;
                    }
                }
            }

            calcVecSystem(n, sumStep, s);
        }
    }

    public static void calcVecSystem(int n, double[] sumStep, Matrix[] s) {
        Matrix A[] = new Matrix[n];

        for (int i = 0 ; i < n ; i++) {
            Matrix sum = new Matrix(n, 1);

            for (int k = i ; k < n ; k++) {

                sum = sum.plus(s[k].mult(sumStep[k]));
            }
            A[i] = sum;
        }

        s[0] = A[0].mult(1 / norm(A[0], n));
        for (int i = 1 ; i < n ; i++) {
            if (sumStep[i] == 0) {
                s[i] = s[i - 1];
            } else {
                Matrix op1 = A[i].mult(pow(norm(A[i - 1], n), 2));
                Matrix op2 = A[i - 1].mult(pow(norm(A[i], n), 2));
                Matrix numerator = op1.minus(op2);
                double denominator = norm(A[i - 1], n) * norm(A[i], n) * sqrt(pow(norm(A[i - 1], n), 2) - pow(norm(A[i], n), 2));
                denominator = 1 / denominator;
                s[i] = numerator.mult(denominator);
            }
        }
    }

    public static Matrix MSH(double epsilon, Matrix x, double m, double b, int n) {
        Matrix xk = x;
        Matrix xk1 = null;

        int counter = 0;
        while (true) {
            System.out.println(counter);
            xk1 = MR1(n, xk, 0.001, m);
            System.out.println("x = [ " + xk1.getByIndex(0, 0) + " , " + xk1.getByIndex(1, 0) + " , " + xk1.getByIndex(2, 0) + "]");
            System.out.println("alpha(x) = " + alpha(xk1));
            System.out.println("m = " + m);


            if (m * alpha(xk1) < epsilon) {
                //System.out.println("x = [ " + xk1.getByIndex(0, 0) + " , " + xk1.getByIndex(1, 0) + "]");

                return xk1;
            }

            m = b * m;
            xk = xk1;

            System.out.println("Falpha(xk) = " + Falpha(xk, m));
            System.out.println("F(xk) = " + fun(xk));
            counter++;
        }
    }

    //matrix as vector
    public static double norm(Matrix x, int n) {
        double sum = 0;
        for (int i = 0 ; i < n ; i++)
            sum += pow(x.getByIndex(i, 0), 2);
        return sqrt(sum);
    }

    public static void main(String[] args) {
        int n = 3;
        Matrix x0 = new Matrix(n, 1);
        x0.setByIndex(0, 0, 7);
        x0.setByIndex(1, 0, 1);
        x0.setByIndex(2, 0, 8);
        //MR1(n, x0, 0.1);

        double m = 0.1;
        double b = 15;
        //MSHGrad(0.1, x0, m, b, n, 0.1, 0.01);
        MSH(0.001, x0, m, b, n);
    }
}
