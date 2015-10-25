

module Science.QuantumChemistry.BasisSet.SerializeBasis where

-- ==================> Standard and third party libraries <=================
import Control.Concurrent.Async
import Control.Monad ((>=>), zipWithM_)
import qualified Data.ByteString as B
import Data.Serialize
import System.Directory (getDirectoryContents)

-- =================> Internal Modules <======================
import Science.QuantumChemistry.ParsecTools.ParserBasis (parseBasisFile)

-- =======================> <==============================

-- | Parse the basis set coefficients and exponents, storing them as a ADT
-- | that is subsequently serialize and store in disk as Bytestring
serializeBasisFile :: FilePath -> IO ()
serializeBasisFile path = do 
   names <- getDirectoryContents path
   let properNames = filter (`notElem` [".", ".."]) names
       basisName   = map ((++ ".basis") . fst . break (== '.') ) properNames
   basis <- mapConcurrently ( parseBasisFile >=> (return . encode )) properNames 
   zipWithM_ B.writeFile basisName basis 
