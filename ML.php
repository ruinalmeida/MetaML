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
	
	function read_data2fromcsv($csv){
		// Read data from data.csv
		$column_x = $column_y = array();
		$f = fopen($csv, 'r');
		while(($data = fgetcsv($f)) !== FALSE){
			$column_x[] = $data[0];
			$column_y[] = $data[1];
		} 
		fclose($f);
		return array("xs"=>$column_x, "ys"=>$column_y);
	}
	
	function gradiant_descent1($data,$learning_rate,$init_bias,$init_slope, $iterations ){	
		$_x=$data['xs'];
		$_y=$data['ys'];
		$error = $this->get_error_for_gd1($bias, $slope, $_x, $_y);
		$bx = $init_bias;
		$mx = $init_slope;
		$num_data_rows = count($_x);

		for($i = 0; $i < $iterations; $i++){
			$_b = 0;
			$_m = 0;
			for($j = 0; $j < $num_data_rows; $j++){
				$x_value = $_x[$j];
				$y_vlaue = $_y[$j];
				$_b += -(2/$num_data_rows) * ($y_vlaue - (($mx * $x_value) + $bx));
				$_m += -(2/$num_data_rows) * ($y_vlaue - (($mx * $x_value) + $bx)) * $x_value;
			}
			$bx = $bx - ($learning_rate * $_b);
			$mx = $mx - ($learning_rate * $_m);
		}
		$final_bias = $bx;
		$final_slope = $mx;
		$error = $this->get_error_for_gd1($current_bias, $current_slope, $_x, $_y);
		return array("bias"=>$final_bias, "slope"=>$final_slope, "error"=>$error,"cx"=>$_x,"cy"=>$_y);
	} 

	function get_error_for_gd1($b, $m, $column_x, $column_y){
		$totalError = 0;
		$num_data_rows = count($column_x);

		for($i = 0; $i < $num_data_rows; $i++){
			$x = $column_x[$i];
			$y = $column_y[$i];
			$totalError += ($y - ($m * $x + $b))*($y - ($m * $x + $b));
		}
		return $totalError / $num_data_rows;
	}
	
}

