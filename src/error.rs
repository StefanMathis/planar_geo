/*!
This module contains the [`Error`] struct, which contains a boxed [`ErrorType`].
The [`ErrorType`] represents all the possible failures which can occur when
using this crate. See the docstrings of struct and enum for details.
*/

use crate::contour::Contour;
use compare_variables::*;
use core::result;

/// An alias for `Result<T, Error>`.
pub type Result<T> = result::Result<T, Error>;

/**
An enum representing all errors which can occur when constructing / modifying a
[`Shape`](crate::shape::Shape).
 */
#[derive(Clone, Debug, PartialEq)]
pub enum ShapeConstructorError<T> {
    /// Contour vector used in [`Shape::new`](crate::shape::Shape::new) is
    /// empty.
    EmptyVec,
    /**
    A shape can only consists of non-empty contours. If [`Contour::is_empty`]
    returns `true` for a contour at index `idx`, the function input and the
    index are returned.
     */
    EmptyContour {
        /// Argument used in the failed function call, can be used for further
        /// inspection.
        input: T,
        /// If `input` was a collection (i.e. a [`Vec<Contour>`]), this value
        /// is the index of the found empty contour.
        idx: usize,
    },
    /**
    An hole is outside the outer contour of the [`Shape`](crate::shape::Shape).
    If `input` is a collection, the index `idx` specifies the hole in question
    and the contour is the first element of the collection.
     */
    HoleOutsideContour {
        /// Argument used in the failed function call, can be used for further
        /// inspection.
        input: T,
        /// If `input` was a collection (i.e. a [`Vec<Contour>`]), this value
        /// is the index of the hole outside the outer contour.
        idx: usize,
    },
    /**
    An hole is inside one of the holes of the [`Shape`](crate::shape::Shape).
    If `input` is a collection, the hole contour at the index `outer_hole_idx`
    contains the hole contour at the index `inner_hole_idx`.
     */
    HoleInsideHole {
        /// Argument used in the failed function call, can be used for further
        /// inspection.
        input: T,
        /// Index of the outer / containing hole
        outer_hole_idx: usize,
        /// Index of the inner / contained hole
        inner_hole_idx: usize,
    },
    /**
    An intersection was found: Either one of the contours of the input
    intersects itself / another contour or one of the existing contours of the
    shape is intersected by the input.
     */
    Intersection {
        /// Argument used in the failed function call, can be used for further
        /// inspection.
        input: T,
        /// Point and indices of the intersection. If the error was returned by
        /// using a function of an existing shape, the left side is the index
        /// of the shape contours and the right side is the index of the
        /// `input`.
        intersection: crate::composite::Intersection,
    },
}

impl<T> std::fmt::Display for ShapeConstructorError<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ShapeConstructorError::EmptyVec => write!(f, "given contour vector was empty"),
            ShapeConstructorError::EmptyContour { input: _, idx } => {
                write!(f, "contour {} is empty", idx)
            }
            ShapeConstructorError::HoleOutsideContour { input: _, idx } => {
                write!(
                    f,
                    "hole contour {} is outside the outer contour of the shape",
                    idx
                )
            }
            ShapeConstructorError::HoleInsideHole {
                input: _,
                outer_hole_idx,
                inner_hole_idx,
            } => write!(
                f,
                "hole contour {} contains the hole at index {}",
                outer_hole_idx, inner_hole_idx
            ),
            ShapeConstructorError::Intersection {
                input: _,
                intersection,
            } => {
                write!(
                    f,
                    "segment {} of contour {} intersects segment {} of contour {} at point x = {}, y = {}",
                    intersection.left.segment_idx,
                    intersection.left.contour_idx,
                    intersection.right.segment_idx,
                    intersection.right.contour_idx,
                    intersection.point[0],
                    intersection.point[1]
                )
            }
        }
    }
}

impl<T: std::fmt::Debug> std::error::Error for ShapeConstructorError<T> {}

/**
This struct represents all errors which can occur when using this crate. This is
a thin wrapper around a boxed [`ErrorType`]. Its main purpose is to provide a
small container (essentially just a pointer) which can be easily passed around
while avoiding moving the large [`ErrorType`] in memory.

See the documentation of [`ErrorType`] for details regarding the different
failure modes.
 */
#[derive(Debug)]
pub struct Error(pub Box<ErrorType>);

impl Error {
    /**
    Returns a reference to the underlying [`ErrorType`].
     */
    pub fn cause(&self) -> &ErrorType {
        return &*self.0;
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        return self.0.fmt(f);
    }
}

impl std::error::Error for Error {}

impl From<ErrorType> for Error {
    fn from(err: ErrorType) -> Self {
        return Error(Box::new(err));
    }
}

impl From<Comparison<f64>> for Error {
    fn from(err: Comparison<f64>) -> Self {
        return ErrorType::Comparison(err).into();
    }
}

/**
An enum containing all errors which can occur when using this crate.
 */
#[derive(Debug)]
pub enum ErrorType {
    /// Received values outside of expected value ranges
    Comparison(Comparison<f64>),
    /// An arc can only be constructed from three points if these points are not
    /// collinear, i.e. are not on a single straight line.
    Collinear([[f64; 2]; 3]),
    /// Tried to construct a segment from (almost) equal start and end points.
    /// The equality is defined by approximate comparison using
    /// [`approx::ulps_eq`], which is why both points are given.
    PointsIdentical {
        /// Start point of the failed segment construction attempt.
        start: [f64; 2],
        /// Stop / end point of the failed segment construction attempt.
        stop: [f64; 2],
    },
    /// Creating a new shape failed because of the contained
    /// [`ShapeConstructorError].
    NewShape(ShapeConstructorError<Vec<Contour>>),
    /// Adding a new hole to an existing shape failed because of the contained
    /// [`ShapeConstructorError].
    AddHole(ShapeConstructorError<Contour>),
}

impl std::fmt::Display for ErrorType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ErrorType::Comparison(err) => err.fmt(f),
            ErrorType::Collinear(pts) => {
                write!(f, "given points {:?} are collinear.", pts)
            }
            ErrorType::PointsIdentical { start, stop } => {
                write!(
                    f,
                    "start point {:?} and end point {:?} are identical",
                    start, stop
                )
            }
            ErrorType::NewShape(err) => {
                return err.fmt(f);
            }
            ErrorType::AddHole(err) => {
                return err.fmt(f);
            }
        }
    }
}

impl std::error::Error for ErrorType {}

impl From<ShapeConstructorError<Vec<Contour>>> for ErrorType {
    fn from(value: ShapeConstructorError<Vec<Contour>>) -> Self {
        return ErrorType::NewShape(value);
    }
}

impl From<ShapeConstructorError<Contour>> for ErrorType {
    fn from(value: ShapeConstructorError<Contour>) -> Self {
        return ErrorType::AddHole(value);
    }
}

impl From<ShapeConstructorError<Vec<Contour>>> for Error {
    fn from(value: ShapeConstructorError<Vec<Contour>>) -> Self {
        return ErrorType::from(value).into();
    }
}

impl From<ShapeConstructorError<Contour>> for Error {
    fn from(value: ShapeConstructorError<Contour>) -> Self {
        return ErrorType::from(value).into();
    }
}
