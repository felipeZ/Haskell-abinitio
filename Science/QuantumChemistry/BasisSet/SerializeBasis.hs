{-|
Module: Science.QuantumChemistry.BasisSet.SerializeBasis
Description: Store plain-text basis into binary format
Copyright: @2016 Felipe Zapata
-}


module Science.QuantumChemistry.BasisSet.SerializeBasis where

-- ==================> Standard and third party libraries <=================
import Control.Concurrent.Async
import Control.Monad ((>=>), zipWithM_)
import qualified Data.ByteString as B
import Data.Serialize
import System.Directory (getDirectoryContents)
import System.FilePath.Posix ((</>))

-- =================> Internal Modules <======================
import Science.QuantumChemistry.ParsecTools.ParserBasis (parseBasisFile)

-- =======================> <==============================

{- | Parse the basis set coefficients and exponents, storing them as a ADT
that is subsequently serialize and store in disk as Bytestring -}
serializeBasisFile :: FilePath -> IO ()
serializeBasisFile path = do 
   names <- getDirectoryContents path
   let properNames = map (path </>) . filter (`notElem` [".", ".."]) $ names
       basisName   = map ((++ ".basis") . takeWhile (/= '.')) properNames
   print properNames
   basis <- mapConcurrently ( parseBasisFile >=> (return . encode )) properNames 
   zipWithM_ B.writeFile basisName basis 
