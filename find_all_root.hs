import Data.Time ( diffUTCTime, getCurrentTime, NominalDiffTime )
import Data.Maybe
import Data.List

-- calc the value of polynomial at point x
-- using foldl as the polynomial is from highest power
polyval ::[Double] -> Double -> Double
polyval coefficients x = foldl (\result coef -> result * x + coef) 0 coefficients


-- calc the derivtive coefficients
-- use map to change int power array to a general type Num
calcDerivPoly :: [Double] -> [Double]
calcDerivPoly coefficients = zipWith (*) (map fromIntegral [length coefficients - 1, length coefficients - 2..1]) coefficients


-- | Calculate the derivative coefficients of a polynomial and normalize them.
--
-- This function takes a list of polynomial coefficients and calculates the coefficients of the derivative polynomial.
calcDerivPolyNormalized :: [Double] -> [Double]
calcDerivPolyNormalized coefficients =
    let n = length coefficients
        powers = [n-1, n-2 .. 1] -- List of powers in reverse order
        normalizedPowers = map (\i -> fromIntegral i / fromIntegral (n-1)) powers -- Normalize the powers
    in zipWith (*) normalizedPowers coefficients

calculatePolynomialDividedByDerivative :: [Double] -> [Double] -> Double -> Maybe Double
calculatePolynomialDividedByDerivative func deriv x =
    let invX = 1/x
        funcV = x * foldr (\result coef -> result * invX + coef) 0 func
        derivV = foldr (\result coef -> result * invX + coef) 0 deriv
    in if derivV == 0
        then Nothing
        else Just (funcV / derivV)


-- | Remove similar values from a list based on a specified tolerance 'epsilon'.
-- This function recursively processes the input list 'xs' and removes values that are considered similar
-- to the preceding value if the absolute difference between them is less than 'epsilon'. It retains the first
-- occurrence of each distinct value and removes subsequent similar values.
removeSimilar :: [Double] -> Double -> [Double]
removeSimilar [] epsilon = []
removeSimilar (x:xs) epsilon
    | null xs = [x]
    | abs (x - head xs) < epsilon =  removeSimilar xs epsilon
    | otherwise = x : removeSimilar xs epsilon


-- | Find the maximum absolute value in a list of numeric values.
-- This function recursively iterates through the input list 'xs' and compares each element's absolute value
-- to the current maximum value 'maxTemp'. If an element with a greater absolute value is found, it updates 'maxTemp'.
--
maxls :: [Double] -> Double -> Double
maxls [] maxTemp = maxTemp
maxls (x:xs) maxTemp
    | maxTemp < abs x = maxls xs (abs x)
    | otherwise = maxls xs maxTemp


-- | Perform a one-directional search on the x-axis until it finds a point where the sign of the polynomial function changes.
-- This function iteratively increments the input 'x' value by the specified 'step' in the direction specified by 'sign'.
-- It evaluates the polynomial represented by 'coefficients' at each step and checks if the sign changes. When a sign change is detected,
-- the function returns the 'x' value at which the sign change occurred.
boundCheck :: [Double] -> Double -> Double -> Double -> Double
boundCheck coefficients x step sign =
    let newX = x + step
    in if signum (polyval coefficients newX) == sign
        then boundCheck coefficients newX step sign
        else newX


--find the maximum abs x value we can find a root
--start boundCheck in the step direction
searchBounds ::  [Double] -> Double -> Double -> Maybe Double
searchBounds coefficients x step =
    let sign = signum (polyval coefficients x)
        --Using Lagrange's Bound 
        maximum_distance = (maxls (tail coefficients) 0.0 / abs (head coefficients)) * fromIntegral (length coefficients) + 1
        headValue = polyval coefficients (sign * maximum_distance)
    in
        if signum (polyval coefficients (signum step * maximum_distance)) == sign
        then Nothing
        else  Just (boundCheck coefficients x step sign)



-- | Perform the Newton-Raphson root-finding method to approximate a root of a polynomial equation.
--
-- Given the coefficients of a polynomial, its derivative coefficients, an initial guess 'x', a desired tolerance 'epsilon',
-- a maximum allowed divergence 'maxDiv', and the maximum number of iterations 'iterations', this function attempts to find
-- an approximate root of the polynomial equation within the specified parameters.
--
-- If the method converges to a root within the given tolerance and divergence constraints, it returns 'Just' the approximate root.
-- If the method fails to find a root within the maximum allowed iterations or due to other conditions, it returns 'Nothing'.
newtonRaphson ::[Double] -> [Double] -> Double -> Double -> Double -> Integer -> Maybe Double
newtonRaphson coefficients dervtive x epsilon maxDiv iterations
    | iterations == 0 = Nothing
    | abs x > 1 =
        let div = calculatePolynomialDividedByDerivative coefficients dervtive x
            xNew = if fromMaybe 0 div == 0
                    then x 
                    else x - fromMaybe 0 div
        in if abs (xNew - x ) > maxDiv||(fromMaybe 0 div == 0)
            then Nothing
            else if  abs (xNew - x ) < epsilon
                then Just xNew
                else newtonRaphson coefficients dervtive xNew epsilon maxDiv (iterations - 1)
    | otherwise =
        let
           func_value = polyval coefficients x
           dervtive_value = polyval dervtive x
           xNew = x - (func_value / dervtive_value)
        in if dervtive_value == 0 ||  abs (xNew - x ) > maxDiv
            then Nothing
            else if abs (xNew - x ) < epsilon
                then Just xNew
                else newtonRaphson coefficients dervtive xNew epsilon maxDiv (iterations - 1)


-- | Perform the bisection method to approximate a root of a polynomial equation.
--
-- Given the coefficients of a polynomial, an initial interval '[a, b]', and a desired tolerance 'epsilon',
-- this function attempts to find an approximate root of the polynomial equation within the specified parameters.
--
-- If the method converges to a root within the given tolerance, it returns the approximate root.
bisection :: [Double] -> Double -> Double -> Double -> Double
bisection coefficients a b epsilon
    | abs (b-a) < epsilon = (b + a) / 2
    | otherwise =
        let mid = (b + a) / 2
            mid_sign = signum (polyval coefficients mid)
            a_sign =  signum (polyval coefficients a)
        in if mid_sign * a_sign > 0
            then bisection coefficients mid b epsilon
            else bisection coefficients a mid epsilon


--find root a in [a,b] using newtonRaphson or bisection
rootFinder :: [Double] -> Double -> Double -> Double ->[Double]
rootFinder coefficients a b epsilon =
    if signum (polyval coefficients a) /= signum (polyval coefficients b)
        then let
            dervtive = calcDerivPoly coefficients
            root = fromMaybe a (newtonRaphson coefficients dervtive ((a+b)/2)  epsilon (abs(b-a)/2) 4)
            in if a < root && root < b
                then [root]
                else let root = bisection coefficients a b epsilon
                     in [root]
    else []


-- findPolynomRoots: finds all polynom roots
findPolynomRoots :: [Double] -> Double -> [Double]
findPolynomRoots coefficients epsilon
  | length coefficients == 2 = [-(head (tail coefficients) / head coefficients)]
  | otherwise =
    let roots = findPolynomRoots (calcDerivPolyNormalized coefficients) epsilon
        leftBound = searchBounds coefficients (head roots) (-0.1)
        rightBound = searchBounds coefficients (last roots) 0.1
        roots' = maybeToList leftBound ++ roots ++ maybeToList rightBound
        result_roots = concatMap
            (\(a, b) -> rootFinder coefficients a b  epsilon)
            (zip roots' (tail roots'))
     in removeSimilar result_roots 0.000001
                                    

measureElapsedTime :: IO a -> IO (a, NominalDiffTime)
measureElapsedTime action = do
  startTime <- getCurrentTime
  result <- action
  endTime <- getCurrentTime
  let elapsedTime = diffUTCTime endTime startTime
  return (result, elapsedTime)


main :: IO ()
main = do
    let epsilon = 1e-12 -- Set the desired precision or tolerance for the roots
    contents <- readFile "poly_coeff(997).txt" -- Read the contents of the file "poly_coeff(997).txt"
    let poly = map read (lines contents) :: [Double] -- Convert the lines to a list of Double coefficients

    -- Measure the elapsed time to find the roots
    (roots, elapsedTime) <- measureElapsedTime ( do
        let roots = findPolynomRoots poly epsilon -- Find the roots of the polynomial using findPolynomRoots
        putStrLn ("The roots from our implementation are:\n" ++ show roots))

    putStrLn ("The elapsed time is: " ++ show elapsedTime) -- Print the elapsed time
    
    -- Evaluate the polynomial at a specific value and print the result
    --let value = -0.2856465079737269
        --value = -0.8315492002712028
        --result = polyval poly value
    --putStrLn ("Result F(-0.2856465079737269) =  " ++ show result)