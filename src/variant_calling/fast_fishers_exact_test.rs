use float_extras::f64;
use std::cmp::{max, min};

const MAXN: i64 = 319795061410741;

fn lgamma(x: i64) -> f64 {
    f64::lgamma(x as f64)
}

fn exp(x: f64) -> f64 {
    x.exp()
}

fn round(x: f64) -> f64 {
    x.round()
}

fn mln_test2t(a: i64, ab: i64, ac: i64, abcd: i64) -> f64 {
    if 0 > a || a > ab || a > ac || ab + ac > abcd + a {
        panic!("invalid contingency table");
    }
    if abcd > MAXN {
        panic!("the grand total of contingency table is too large");
    }
    let a_min = max(0, ab + ac - abcd);
    let a_max = min(ab, ac);
    if a_min == a_max {
        return 0.0;
    }
    let p0 = lgamma(ab + 1) + lgamma(ac + 1) + lgamma(abcd - ac + 1) + lgamma(abcd - ab + 1)
        - lgamma(abcd + 1);
    let pa =
        lgamma(a + 1) + lgamma(ab - a + 1) + lgamma(ac - a + 1) + lgamma(abcd - ab - ac + a + 1);
    let mut st = 1.0;
    if ab * ac < a * abcd {
        for i in (a_min..min(a, round(ab as f64 * ac as f64 / abcd as f64) as i64 + 1)).rev() {
            let pi = lgamma(i + 1)
                + lgamma(ab - i + 1)
                + lgamma(ac - i + 1)
                + lgamma(abcd - ab - ac + i + 1);
            if pi < pa {
                continue;
            }
            let st_new = st + exp(pa - pi);
            if st_new == st {
                break;
            }
            st = st_new;
        }
        for i in a + 1..a_max + 1 {
            let pi = lgamma(i + 1)
                + lgamma(ab - i + 1)
                + lgamma(ac - i + 1)
                + lgamma(abcd - ab - ac + i + 1);
            let st_new = st + exp(pa - pi);
            if st_new == st {
                break;
            }
            st = st_new;
        }
    } else {
        for i in (a_min..a).rev() {
            let pi = lgamma(i + 1)
                + lgamma(ab - i + 1)
                + lgamma(ac - i + 1)
                + lgamma(abcd - ab - ac + i + 1);
            let st_new = st + exp(pa - pi);
            if st_new == st {
                break;
            }
            st = st_new;
        }
        for i in max(a + 1, round(ab as f64 * ac as f64 / abcd as f64) as i64)..a_max + 1 {
            let pi = lgamma(i + 1)
                + lgamma(ab - i + 1)
                + lgamma(ac - i + 1)
                + lgamma(abcd - ab - ac + i + 1);
            if pi < pa {
                continue;
            }
            let st_new = st + exp(pa - pi);
            if st_new == st {
                break;
            }
            st = st_new;
        }
    }
    (pa - p0 - st.ln()).max(0.0)
}

pub fn fishers_exact_test(a: i64, b: i64, c: i64, d: i64) -> f64 {
    exp(-mln_test2t(a, a + b, a + c, a + b + c + d))
}
