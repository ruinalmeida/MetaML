<?php
namespace metaml;
class Matrix {
	
	function matrix_multiply($m1,$m2){
		$r=count($m1);
		$c=count($m2[0]);
		$p=count($m2);
		if(count($m1[0])!=$p){throw new \Exception('Invalid dimensions');}
		$m3=array();
		for ($i=0;$i< $r;$i++){
			for($j=0;$j<$c;$j++){
				$m3[$i][$j]=0;
				for($k=0;$k<$p;$k++){
					$m3[$i][$j]+=$m1[$i][$k]*$m2[$k][$j];
				}
			}
		}
		return($m3);
	}

	function matrix_sum($m1,$m2){
		$r=count($m1);
		$c=count($m2[0]);
		$p=count($m2);
		$t=count($m1[0]);
		if($r!=$p || $c!=$t ){throw new \Exception('Invalid dimensions');}
		$m3=array();
		for ($i=0;$i< $r;$i++){
			for($j=0;$j<$c;$j++){
				$m3[$i][$j]=$m1[$i][$j] +$m2[$i][$j];
			}
		}
		return($m3);
	}
	
	function matrix_subtract($m1,$m2){
		$r=count($m1);
		$c=count($m2[0]);
		$p=count($m2);
		$t=count($m1[0]);
		if($r!=$p || $c!=$t ){throw new \Exception('Invalid dimensions');}
		$m3=array();
		for ($i=0;$i< $r;$i++){
			for($j=0;$j<$c;$j++){
				$m3[$i][$j]=$m1[$i][$j] -$m2[$i][$j];
			}
		}
		return($m3);
	}
	
	function get_identity($n){
		$m3=array();
		for ($i=0;$i< $n;$i++){
			for($j=0;$j<$n;$j++){
				if ($i==$j)
				{
					$m3[$i][$j]=1;	
				}
				else
				{
					$m3[$i][$j]=0;
				}

			}
		}
		return($m3);
	}
	
	function matrix_transppose($m){
		$r=count($m);
		$c=count($m[0]);
		$mt=array();
		for($i=0;$i< $r;$i++){
			for($j=0;$j<$c;$j++){
				$mt[$j][$i]=$m[$i][$j];
			}
		}
		return($mt);
	}
	
	function matrix_determinant($m) {
			$r=count($m);
			$c=count($m[0]);
		    $norm = 0;
			$count=0;
			if ($r != $c) {
					throw new \Exception('Invalid m'.$c.'dimensions'.$r);
            }
			// faz a triangular transform...
            for ($i=0; $i < $r-1; $i++) {

				if($m[$i][$i] == 0)
				{
					for($k = $i; $k < $r; $k++)
					{
						if($m[$k][$i] != 0)
						{
							for($j = 0; $j < $r; $j++)
							{
								$temp = $m[$i][$j];
								$m[$i][$j] = $m[$k][$j];
								$m[$k][$j] = $temp;
							}
							$k = $r;
						}
					}
					$count++;
				}
				if($m[$i][$i] != 0)
					{
						for($k = $i + 1; $k < $r; $k++)
						{
							$factor = -1.0 * $m[$k][$i] /  $m[$i][$i];
							for($j = $i; $j < $r; $j++)
							{
								$m[$k][$j] = $m[$k][$j] + ($factor * $m[$i][$j]);
							}
						}
					}
            }
			$temp = 1.0;
			// get determinant
			for($i = 0; $i < $r; $i++)
			{
				$temp = $temp*$m[$i][$i];
			}
           if($count % 2 == 0)
		   {
			   return $temp;
		   }
        
			else
			{
				return   -1.0 * $temp;
			}

	}

	// LU factorization with pivoting.
    // perm are row permutations; toggle is +1 or -1 (even or odd)
	function matrix_decompose($matrix, &$perm,&$toggle){
			$rows=count($matrix);
			$cols=count($matrix[0]);
            if ($rows != $cols)
			{
				throw new \Exception('Invalid dimensions ');
			}
                
            $n = $rows; 

			$perm = array();
            for ($i = 0; $i < $n; ++$i) 
			{ 
				$perm[$i] = $i; 
			}

            $toggle = 1; 
			
			for ($j = 0; $j < $n - 1; ++$j) 
            {
                $colMax = sqrt($matrix[$j][$j]*$matrix[$j][$j]); 
                $pRow = $j;

                for ($i = $j + 1; $i < $n; ++$i)
                {
                    if (sqrt($matrix[$i][$j]*$matrix[$i][$j]) > $colMax)
                    {
                        $colMax = sqrt($matrix[$i][$j]*$matrix[$i][$j]);
                        $pRow = $i;
                    }
                }
         
                if ($pRow != $j) 
                {
                    $rowPtr = $matrix[$pRow];
                    $matrix[$pRow] = $matrix[$j];
                    $matrix[$j] = $rowPtr;

                    $tmp = $perm[$pRow]; 
                    $perm[$pRow] = $perm[$j];
                    $perm[$j] = $tmp;

                    $toggle = -$toggle; 
                }

                if ($matrix[$j][$j] == 0.0)
                {
                 
                    $goodRow = -1;
                    for ($row = $j + 1; $row < $n; ++$row)
                    {
                        if ($matrix[$row][$j] != 0.0)
                            $goodRow = $row;
                    }

                    if ($goodRow == -1)
                        new \Exception("Cannot use Doolittle's method");

                   
                    $rowPtr = $matrix[$goodRow];
                    $matrix[$goodRow] = $matrix[$j];
                    $matrix[$j] = $rowPtr;

                    $tmp = $perm[$goodRow]; 
                    $perm[$goodRow] = $perm[$j];
                    $perm[$j] = $tmp;

                    $toggle = -$toggle; 
                }
  
                for ($i = $j + 1; $i < $n; ++$i)
                {
                    $matrix[$i][$j] /= $matrix[$j][$j];
                    for ($k = $j + 1; $k < $n; ++$k)
                    {
                        $matrix[$i][$k] -= $matrix[$i][$j] * $matrix[$j][$k];
                    }
                }

            } 
			return $matrix;
			
		}
	
	function helper_solve($lu, $h){
            $n = count($lu);
            for ($i = 1; $i < $n; ++$i)
            {
				$sum = $h[$i];
                for ($j = 0; $j < $i; ++$j)
				{
					$sum -= $lu[$i][$j] * $h[$j];
					$h[$i] = $sum;
				}
                    
            }
            $h[$n - 1] /= $lu[$n - 1][$n - 1];
            for ($i = $n - 2; $i >= 0; --$i)
            {
                $sum = $h[$i];
                for ($j = $i + 1; $j < $n; ++$j)
				{
					$sum -= $lu[$i][$j] * $h[$j];
					$h[$i] = $sum / $lu[$i][$i];	
				}

            }
            return $h;
        }
    
	function matrix_inverse($matrix){
            $n = count($matrix);

            $lu = $this->matrix_decompose($matrix, $perm,$toggle);
            if ($lu == null)
                throw new \Exception("Unable to compute inverse");

            $h = array();
			//  permute h using the permutation array used in LU decompose
            for ($i = 0; $i < $n; ++$i)
            {
                for ($j = 0; $j < $n; ++$j)
                {
                    if ($i == $perm[$j])
					{
						$h[$j] = 1.0;
					} 
                    else
					{
						$h[$j] = 0.0;
					}    
                }
                $x = $this->helper_solve($lu, $h); 
                for ($j = 0; $j < $n; ++$j)
				{
					$matrix[$j][$i] = $x[$j];
				}      
            }
            return $matrix;
        }

}
