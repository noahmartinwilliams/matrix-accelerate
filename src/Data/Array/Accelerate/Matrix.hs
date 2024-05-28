{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
{-|
Module        : Data.Array.Accelerate.Matrix
Description   : Functions for plain and dependently typed matrix math.
Copyright     : (c) Noah Martin Williams 2024
License       : BSD-3-Clause
Maintainer    : noahmartinwilliams@gmail.com
Stability     : experimental
Portability   : POSIX

This module contains functions for doing matrix math such as addition, subtraction, and multiplication
for both plain and dependently typed matrices.
-}
module Data.Array.Accelerate.Matrix(mMul, matMul, identMat, Mat(..), AccMat(..), matTransp, matAdd, mAdd, mSub, matSub, useMat) where

import Prelude as P
import Data.Array.Accelerate as A
import Data.Array.Accelerate.Control.Lens.Shape

-- |Creat an identity matrix with the dimension provided
identMat :: Exp Int -> Acc (Matrix Int)
identMat n = generate (index2 n n) (\(I2 a b) -> a A.== b A.? ((constant 1), (constant 0)))

-- These two lines were blatantly stolen from here: https://hackage.haskell.org/package/accelerate-1.3.0.0/docs/Data-Array-Accelerate.html
rep0 :: (Shape sh, Elt e) => Exp Int -> Acc (Array sh e) -> Acc (Array (sh:.Int) e)
rep0 n a = A.replicate (lift (Any:.n)) a

{-
    Example:
    let m1 = A.generate (constant (Z:.3:.5)) (\(I2 x y) -> (A.fromIntegral x :: Exp Double) * 5.0 + (A.fromIntegral y :: Exp Double))
        m2 = A.generate (constant (Z:.5:.3)) (\(I2 x y) -> (A.fromIntegral x :: Exp Double) * 5.0 + (A.fromIntegral y :: Exp Double))
        m1 `mMul` m2 = 
   [150   160   170
    400   435   470
    650   710   770]
-}

-- |Multiply two matrices together without dependent types.
mMul :: A.Num e => Acc (Matrix e) -> Acc (Matrix e) -> Acc (Matrix e)
mMul m1 m2 = do
    let Z:.m1NumRows:._ = unlift (shape m1 :: Exp (Plain (Z:.Int:.Int))) :: (Z:.Exp Int:.Exp Int)
        Z:._:.m2NumCols = unlift (shape m2 :: Exp (Plain (Z:.Int:.Int))) :: (Z:.Exp Int:.Exp Int)
        m1Tensor = transposeOn _1 _2 (rep0 m2NumCols m1)
        m2Tensor = transposeOn _2 _3 (transposeOn _1 _2 (rep0 m1NumRows (A.transpose m2)))
        ret = A.sum (A.zipWith (*) m1Tensor m2Tensor)
    ret

-- |Add two matrices without dependent types.
mAdd :: A.Num e => Acc (Matrix e) -> Acc (Matrix e) -> Acc (Matrix e)
mAdd left right = A.zipWith (+) left right

-- |Subtract one matrix from another without dependent types.
mSub :: A.Num e => Acc (Matrix e) -> Acc (Matrix e) -> Acc (Matrix e)
mSub left right = A.zipWith (-) left right

-- |Dependent type for accelerate matrices.
data AccMat e a b where 
    -- |Dependently typed accelerated matrix which forces two types to line up.
    AccMat :: (Elt e, A.Num e) => Acc (Matrix e) -> a -> b -> AccMat e a b

-- |Dependent type for plain matrices.
data Mat e a b where
    -- |Dependently typed plain matrix for passing to compiled functions.
    Mat :: (Elt e, A.Num e) => Matrix e -> a -> b -> Mat e a b

-- |Change the type of a dependently typed matrix from AccMat to Mat.
useMat :: Mat e a b -> AccMat e a b
useMat (Mat mat a b) = AccMat (use mat) a b

-- |Multiply two dependently typed matrices together.
-- |For example:
-- 
-- @
-- data A = A
-- data B = B
-- data C = C
-- 
-- let m1 = AccMat (use (fromList (Z:.10:.12) [0..] :: Matrix Int)) A B
-- let m2 = AccMat (use (fromList (Z:.12:.13) [0..] :: Matrix Int)) B C
-- let mResult = m1 `matMul` m2
matMul :: A.Num e => AccMat e a b -> AccMat e b c -> AccMat e a c
matMul (AccMat left a _) (AccMat right _ c) = AccMat (mMul left right) a c

-- |Add two dependently typed matrices.
matAdd :: A.Num e => AccMat e a b -> AccMat e a b -> AccMat e a b
matAdd (AccMat left a b) (AccMat right _ _) = AccMat (mAdd left right) a b

-- |Subtract one dependently typed matrix from another.
matSub :: A.Num e => AccMat e a b -> AccMat e a b -> AccMat e a b
matSub (AccMat left a b) (AccMat right _ _) = AccMat (mSub left right) a b

-- |Transpose a dependently typed matrix.
matTransp :: AccMat e a b -> AccMat e b a
matTransp (AccMat mat a b) = AccMat (A.transpose mat) b a
