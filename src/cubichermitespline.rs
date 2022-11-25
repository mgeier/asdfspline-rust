use crate::utilities::GridError;
use crate::PiecewiseCubicCurve;
use crate::Vector;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("There must be at least two positions")]
    LessThanTwoPositions,
    #[error(
        "Exactly 2 tangents per segment are required \
            (got {segments} segments and {tangents} tangents)"
    )]
    TangentsVsSegments { tangents: usize, segments: usize },
    #[error("length of grid ({grid}) must be the same as number of positions ({positions})")]
    GridVsPositions { grid: usize, positions: usize },
    #[error(transparent)]
    FromGridError(#[from] GridError),
}

impl<V: Vector> PiecewiseCubicCurve<V> {
    pub fn new_hermite(
        positions: &[V],
        tangents: &[V],
        grid: &[f32],
    ) -> Result<PiecewiseCubicCurve<V>, Error> {
        use Error::*;
        if positions.len() < 2 {
            return Err(LessThanTwoPositions);
        }
        let segments_len = positions.len() - 1;
        if tangents.len() != 2 * segments_len {
            return Err(TangentsVsSegments {
                tangents: tangents.len(),
                segments: segments_len,
            });
        }
        if positions.len() != grid.len() {
            return Err(GridVsPositions {
                grid: grid.len(),
                positions: positions.len(),
            });
        }
        let mut segments = Vec::with_capacity(segments_len);
        for i in 0..segments_len {
            let t0 = grid[i];
            let t1 = grid[i + 1];
            let delta = t1 - t0;
            let x0_delta = positions[i] / delta;
            let x1_delta = positions[i + 1] / delta;
            let v0 = tangents[2 * i];
            let v1 = tangents[2 * i + 1];

            // [a0]                     [x0 / delta]
            // [a1] = delta**(-2) * M * [x1 / delta]
            // [a2]                     [v0        ]
            // [a3]                     [v1        ]

            // M =
            //
            // [t1**2*(-3*t0 + t1), t0**2*(-t0 + 3*t1),      -t0*t1**2,      -t0**2*t1]
            // [           6*t0*t1,           -6*t0*t1, t1*(2*t0 + t1), t0*(t0 + 2*t1)]
            // [      -3*t0 - 3*t1,        3*t0 + 3*t1,     -t0 - 2*t1,     -2*t0 - t1]
            // [                 2,                 -2,              1,              1]

            let t0_2 = t0.powi(2);
            let t1_2 = t1.powi(2);

            segments.push([
                          delta.powi(-2) * (
            x0_delta * t1_2 * (-3*t0 + t1) + x1_delta * t0_2*(-t0 + 3*t1) + v0 * -t0*t1_2 + v1 *   -t0_2*t1
            ),
                          delta.powi(-2) * (
             x0_delta *          6*t0*t1 + x1_delta *           -6*t0*t1 + v0 *  t1*(2*t0 + t1) + v1 *  t0*(t0 + 2*t1)
            ),
                          delta.powi(-2) * (
              x0_delta *    -3*t0 - 3*t1 + x1_delta *        3*t0 + 3*t1 + v0 *      -t0 - 2*t1 + v1 *      -2*t0 - t1
            ),
                          delta.powi(-2) * (
                     x0_delta *        2 + x1_delta *                 -2 + v0 *               1 + v1 *               1
            ),
            ]);
        }
        use crate::piecewisecubiccurve::Error as Other;
        PiecewiseCubicCurve::new(segments, grid).map_err(|err| match err {
            Other::ZeroSegments => unreachable!(),
            Other::GridVsSegments { .. } => unreachable!(),
            Other::FromGridError(e) => e.into(),
        })
    }
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    use crate::Spline; // for grid()

    #[test]
    fn test_1d() {
        let positions = [1.0, 2.0];
        let tangents = [3.0, 4.0];
        let grid = [5.0, 6.0];
        let curve = PiecewiseCubicCurve::new_hermite(&positions, &tangents, &grid).unwrap();
        assert_eq!(curve.grid(), &[5.0, 6.0]);
    }

    #[test]
    fn test_errors() {
        let result = PiecewiseCubicCurve::new_hermite(&[1.0], &[3.0, 4.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = PiecewiseCubicCurve::new_hermite(&[1.0, 2.0], &[3.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = PiecewiseCubicCurve::new_hermite(&[1.0, 2.0], &[3.0, 3.5, 4.0], &[5.0, 6.0]);
        assert!(result.is_err());
        let result = PiecewiseCubicCurve::new_hermite(&[1.0, 2.0], &[3.0, 4.0], &[5.0]);
        assert!(result.is_err());
    }
}
