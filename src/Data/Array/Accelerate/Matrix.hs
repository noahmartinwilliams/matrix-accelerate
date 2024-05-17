{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs #-}
module Data.Array.Accelerate.Matrix(mMul, matMul, identMat, Mat(..), AccMat(..)) where

import Prelude as P
import Data.Array.Accelerate as A
import Data.Array.Accelerate.LLVM.PTX
import Data.Array.Accelerate.Control.Lens.Shape

identMat :: Exp Int -> Acc (Matrix Int)
identMat n = generate (index2 n n) (\(I2 a b) -> a A.== b A.? ((constant 1), (constant 0)))

-- These two lines were blatantly stolen from here: https://hackage.haskell.org/package/accelerate-1.3.0.0/docs/Data-Array-Accelerate.html
rep0 :: (Shape sh, Elt e) => Exp Int -> Acc (Array sh e) -> Acc (Array (sh:.Int) e)
rep0 n a = A.replicate (lift (Any:.n)) a

{-
    Example:
    let m1 = A.generate (constant (Z:.3:.5)) (\(I2 x y) -> (A.fromIntegral x :: Exp Double) * 5.0 + (A.fromIntegral y :: Exp Double))
        m2 = A.generate (constant (Z:.5:.3)) (\(I2 x y) -> (A.fromIntegral x :: Exp Double) * 5.0 + (A.fromIntegral y :: Exp Double))
        m1 * m2 = 
   [150   160   170
    400   435   470
    650   710   770]
-}

mMul :: (Elt e, A.Num e) => Acc (Matrix e) -> Acc (Matrix e) -> Acc (Matrix e)
mMul m1 m2 = do
    let Z:.m1NumRows:._ = unlift (shape m1 :: Exp (Plain (Z:.Int:.Int))) :: (Z:.Exp Int:.Exp Int)
        Z:._:.m2NumCols = unlift (shape m2 :: Exp (Plain (Z:.Int:.Int))) :: (Z:.Exp Int:.Exp Int)
        m1Tensor = transposeOn _1 _2 (rep0 m2NumCols m1)
        m2Tensor = transposeOn _2 _3 (transposeOn _1 _2 (rep0 m1NumRows (A.transpose m2)))
        ret = A.sum (A.zipWith (*) m1Tensor m2Tensor)
    ret

data AccMat e a b where 
    AccMat :: (Elt e, A.Num e) => Acc (Matrix e) -> a -> b -> AccMat e a b

data Mat e a b where
    Mat :: (Elt e, A.Num e) => Matrix e -> a -> b -> Mat e a b

matMul :: (Elt e, A.Num e) => AccMat e a b -> AccMat e b c -> AccMat e a c
matMul (AccMat left a _) (AccMat right _ c) = AccMat (mMul left right) a c

matTransp :: AccMat e a b -> AccMat e b a
matTransp (AccMat mat a b) = AccMat (A.transpose mat) a b
