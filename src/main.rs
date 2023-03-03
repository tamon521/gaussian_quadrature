use std::f64::consts::PI;
fn main() {
    let n = 5;
    let f = |x: f64| x.exp(); // 関数f(x)を定義
    let result = gaussian_quadrature(n, f); // 関数f(x)をn次のGaussian Quadratureで積分
    println!("The result          = {}", result);
    println!("The analysis result = {}", 1.0_f64.exp() - (-1.0_f64).exp());
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
    // let m = (n + 1) / 2;
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    // let eps = 1e-15;
    let mut z;
    let mut x0;
    for i in 0..n {
        x0 = ((i as f64 + 0.75) / (n as f64 + 0.5) * PI).cos();
        z = zeropoint_legendre_poly(n, x0);
        x[i] = z;
        w[i] = 2.0 * (1.0 - z * z)
            / (n as f64 * legendre_poly(n - 1, z))
            / (n as f64 * legendre_poly(n - 1, z));
    }
    (x, w)
}

#[test]
fn test_gauss_legendre() {
    assert_eq!(gauss_legendre(1), (vec![0.0], vec![2.0]));
    assert_eq!(
        gauss_legendre(2),
        (
            vec![1.0_f64 / 3.0_f64.sqrt(), -1.0_f64 / 3.0_f64.sqrt()],
            vec![1.0; 2]
        )
    );
    assert_eq!(
        gauss_legendre(3),
        (
            vec![0.6_f64.sqrt(), 0.0, -0.6_f64.sqrt()],
            vec![5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
        )
    );
    // todo ここの処理を書く + fn gauss_legendre を解読して動くようにする
}

fn zeropoint_legendre_poly(n: usize, mut x0: f64) -> f64 {
    let imax: i32 = 200;
    for i in 1..=imax {
        let x1 = x0 - legendre_poly(n, x0) / diff_legendre_poly(n, x0);
        if (x1 - x0).abs() < 1e-16 {
            break;
        }
        x0 = x1;
        if i == imax {
            println!("zero point not found");
        }
    }
    if x0.abs() < 1e-10 {
        x0 = 0.0;
    }
    x0
}
#[test]
fn test_zeropoint_legendre_poly() {
    assert!((zeropoint_legendre_poly(1, 0.5) - 0.0).abs() < 1e-10);
    assert!((zeropoint_legendre_poly(2, 0.5) - 1.0 / 3.0_f64.sqrt()).abs() < 1e-10);
}

fn diff_legendre_poly(n: usize, mut x: f64) -> f64 {
    // let eps = 1e-10;
    // if (x - 1.0).abs() < eps {
    //     return (n * (n + 1)) as f64 / 2.0;
    // };
    // if (x + 1.0).abs() < eps {
    //     return power(-1.0, (n + 1).into()) * (n * (n + 1)) as f64 / 2.0;
    // };
    if n == 0 {
        return 0.0;
    }
    if (x == 1.0) || (x == -1.0) {
        x += 0.0001_f64;
    }
    return (n as f64) * (legendre_poly(n - 1, x) - x * legendre_poly(n, x)) / (1.0 - x * x);
}
#[test]
fn test_diff_legendre_poly() {
    assert_eq!(diff_legendre_poly(1, 0.0), 1.0);
    assert_eq!(diff_legendre_poly(2, 0.0), 0.0);
    assert_eq!(diff_legendre_poly(2, 0.5), 1.5);
    assert_eq!(diff_legendre_poly(4, 0.5), -1.5625);
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
#[test]
fn test_legendre_poly() {
    assert_eq!(legendre_poly(1, 0.0), 0.0);
    assert_eq!(legendre_poly(2, 0.0), -0.5);
    assert_eq!(legendre_poly(2, 2.0), 5.5);
    assert_eq!(legendre_poly(3, 1.0), 1.0);
    assert_eq!(legendre_poly(10, 1.0), 1.0);
}

// fn power(n: f64, m: i16) -> f64 {
//     match m {
//         0 => 1.0,
//         1 => n,
//         _ => {
//             let x = power(n, m / 2);
//             if m % 2 == 0 {
//                 x * x
//             } else {
//                 x * x * n
//             }
//         }
//     }
// }
