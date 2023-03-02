use std::f64::consts::PI;
fn main() {
    let n = 10;
    let f = |x: f64| x.exp(); // 関数f(x)を定義
    let result = gaussian_quadrature(n, f); // 関数f(x)をn次のGaussian Quadratureで積分
    println!("The result is {}", result);
}

fn gaussian_quadrature<F>(n: usize, f: F) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut sum = 0.0;
    let (x, w) = gauss_legendre(n);
    for i in 0..n {
        sum += w[i] * f(x[i]);
    }
    sum
}

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    let m = (n + 1) / 2;
    for i in 0..m {
        let z = ((i as f64) + 0.75) * PI / (n as f64 + 0.5);
        let mut p = 1.0;
        let mut pp = 0.0;
        for _j in 0..n {
            let tmp = p;
            p = legendre_poly(_j, z) * (2.0 / (n as f64) + 1.0) - pp / (n as f64) + 1.0;
            pp = tmp;
        }
        x[i] = -p.sqrt();
        x[n - 1 - i] = -x[i];
        w[i] = 2.0 / ((1.0 - x[i] * x[i]) * pp * pp);
        w[n - 1 - i] = w[i];
    }
    (x, w)
}

fn legendre_poly(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut p0 = 1.0;
    let mut p1 = x;
    for i in 1..n {
        let p = ((2.0 * i as f64 + 1.0) * x * p1 - i as f64 * p0) / (i as f64 + 1.0);
        p0 = p1;
        p1 = p;
    }
    p1
}
