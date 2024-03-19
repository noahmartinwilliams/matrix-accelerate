{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
module Data.Array.Accelerate.Matrix where

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
    let m1 = A.generate (constant (Z:.3:.5)) (\(I2 x y) -> (A.fromIntegral x :: Exp Double) * 5.0 + (A.fromIntegral y :: Exp Double))
        m2 = A.generate (constant (Z:.5:.3)) (\(I2 x y) -> (A.fromIntegral x :: Exp Double) * 5.0 + (A.fromIntegral y :: Exp Double))
        m1 * m2 = 
   [150   160   170
    400   435   470
    650   710   770]
-}
matMul :: (Elt e, A.Num e) => Acc (Matrix e) -> Acc (Matrix e) -> Acc (Matrix e)
matMul m1 m2 = do
    let Z:.m1NumRows:._ = unlift (shape m1 :: Exp (Plain (Z:.Int:.Int))) :: (Z:.Exp Int:.Exp Int)
        Z:._:.m2NumCols = unlift (shape m2 :: Exp (Plain (Z:.Int:.Int))) :: (Z:.Exp Int:.Exp Int)
        m1Tensor = transposeOn _1 _2 (rep0 m2NumCols m1)
        m2Tensor = transposeOn _2 _3 (transposeOn _1 _2 (rep0 m1NumRows (A.transpose m2)))
        ret = A.sum (A.zipWith (*) m1Tensor m2Tensor)
    ret
