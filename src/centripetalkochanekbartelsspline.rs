use num_traits::{one, pow, zero};

use crate::make_cubic_hermite_spline;
use crate::PiecewiseCubicCurve;
use crate::{fail, Error};
use crate::{Scalar, Vector};

// TODO: require alga::linear::NormedSpace to get norm() method for vectors?
pub fn make_centripetal_kochanek_bartels_spline<S, V, F>(
    positions: &[V],
    tcb: &[[S; 3]],
    closed: bool,
    get_length: F,
) -> Result<PiecewiseCubicCurve<S, V>, Error>
where
    S: Scalar,
    V: Vector<S>,
    F: Fn(V) -> S,
{
    if positions.len() < 2 {
        fail!("At least two positions are required");
    }
    let mut positions = positions;
    // Only used for "closed" splines:
    let mut positions_vec;
    if closed {
        positions_vec = Vec::with_capacity(positions.len() + 2);
        positions_vec.extend(positions);
        positions_vec.push(positions[0]);
        positions_vec.push(positions[1]);
        positions = &positions_vec;
    } else {
        // To avoid error: "use of possibly uninitialized `positions_vec`"
        positions_vec = Vec::new();
    }

    if tcb.len() + 2 != positions.len() {
        fail!("There must be two more positions than TCB values (except for closed curves)");
    }

    // Create grid with centripetal parametrization

    let mut grid = Vec::<S>::with_capacity(positions.len());
    grid.push(zero());
    for i in 0..positions.len() - 1 {
        let x0 = positions[i];
        let x1 = positions[i + 1];
        let delta = get_length(x1 - x0).sqrt();
        if delta == zero() {
            fail!("Repeated positions are not possible");
        }
        grid.push(*grid.last().unwrap() + delta);
    }
    let mut tangents = Vec::<V>::new();
    assert!(positions.len() == grid.len());
    assert!(positions.len() == tcb.len() + 2);
    for i in 0..positions.len() - 2 {
        let x_1 = positions[i];
        let x0 = positions[i + 1];
        let x1 = positions[i + 2];
        let t_1 = grid[i];
        let t0 = grid[i + 1];
        let t1 = grid[i + 2];
        #[allow(non_snake_case)]
        let [T, C, B] = tcb[(i + closed as usize) % tcb.len()];
        let one = one::<S>();
        let a = (one - T) * (one + C) * (one + B);
        let b = (one - T) * (one - C) * (one - B);
        let c = (one - T) * (one - C) * (one + B);
        let d = (one - T) * (one + C) * (one - B);

        let incoming = ((x0 - x_1) * c * pow(t1 - t0, 2) + (x1 - x0) * d * pow(t0 - t_1, 2))
            / ((t1 - t0) * (t0 - t_1) * (t1 - t_1));
        let outgoing = ((x0 - x_1) * a * pow(t1 - t0, 2) + (x1 - x0) * b * pow(t0 - t_1, 2))
            / ((t1 - t0) * (t0 - t_1) * (t1 - t_1));
        tangents.push(incoming);
        tangents.push(outgoing);
    }

    if closed {
        // Move last (outgoing) tangent to the beginning:
        tangents.rotate_right(1);

        // Remove temporary position and grid elements:
        positions_vec.pop();
        grid.pop();

        // Update reference
        positions = &positions_vec;
    } else if positions.len() == 2 {
        // Straight line
        assert!(grid.len() == 2);
        assert!(tangents.is_empty());
        let tangent = (positions[1] - positions[0]) / (grid[1] - grid[0]);
        tangents.push(tangent);
        tangents.push(tangent);
    } else {
        // End conditions for non-closed curves
        assert!(tangents.len() >= 2);

        let one: S = one();
        let two = one + one;
        let three = two + one;

        // "natural" end conditions
        let natural_end_tangent = |x0, x1, t0, t1, inner_tangent| {
            let delta = t1 - t0;
            (x1 * three - x0 * three - inner_tangent * delta) / (two * delta)
        };

        tangents.insert(
            0,
            natural_end_tangent(positions[0], positions[1], grid[0], grid[1], tangents[0]),
        );
        let last = positions.len() - 1;
        tangents.push(natural_end_tangent(
            positions[last - 1],
            positions[last],
            grid[last - 1],
            grid[last],
            *tangents.last().unwrap(),
        ));
    }
    make_cubic_hermite_spline(&positions, &tangents, &grid)
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn test_1d() {
        let positions = [1.0f32, 2.0, 3.0].to_vec();
        let tcb = [[4.0, 5.0, 6.0]];
        let closed = false;
        let curve = make_centripetal_kochanek_bartels_spline(&positions, &tcb, closed, |x| x.abs())
            .unwrap();
        assert_eq!(curve.grid()[0], 0.0);
        assert_eq!(curve.evaluate(0.0), 1.0);
        assert_eq!(curve.evaluate(*curve.grid().last().unwrap()), 3.0);
    }
}
