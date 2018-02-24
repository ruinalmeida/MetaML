<?php
namespace metaml;

class ML {
	
	/**
	 * sipple linear regression 
	 * @param $xs array x-coords 
	 * @param $ys array y-coords
	 * @returns array() b=>slope, a=>intercept
	 * Y=bX+a
		Regression Equation(y) = a + bx 
		Slope(b) = (NΣXY - (ΣX)(ΣY)) / (NΣX2 - (ΣX)2)
		Intercept(a) = (ΣY - b(ΣX)) / N
	 */
	function linear_regression($xs, $ys) {
		$n = count($xs);
		if ($n != count($ys)) {
			throw new \Exception('Invalid dimensions');
		}
		// calculate sums
		$sumX = array_sum($xs);
		$sumY = array_sum($ys);

		$sumXX = 0;
		$sumXY = 0;
		  
		for($i = 0; $i < $n; $i++) {
    $sumXY+=($xs[$i]*$ys[$i]);
    $sumXX+=($xs[$i]*$xs[$i]);
		}
		//slope
		$b = (($n * $sumXY) - ($sumX * $sumY)) / (($n * $sumXX) - ($sumX * $sumX));
		//intercept
		$a = ($sumY - ($b * $sumX)) / $n;
		// return result
		return array("b"=>$b, "a"=>$a);
	}
	
}
