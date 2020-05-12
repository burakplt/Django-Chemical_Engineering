  function searchBox(){  
    $("select:not(#id_compound_no,#id_kij,#unit-pressure,#unit-temp)").select2( {
     placeholder: "Select",
     allowClear: false
     } )
    }
    
    function convertPressure(value, base, target){
      var factors = new Array(1, 1000, 6894.757, 100000, 101325, 133.3224);
      var result = value;
      result = result * factors[base];
      result = result / factors[target];
      return result;
    }

    function convertTemperature(value, base, target){
        var factors = new Array(1, 1, 0.55555555, 0.5555555)
        var increment = new Array(0, 273.15, 459.67, 0)
        var result = value + increment[base];
        result = result*factors[base];
        result = result/factors[target] - increment[target];
        return result;
    }
  
    searchBox()
    var $unitPressure = $("#unit-pressure"); $basePressure = parseInt($unitPressure.val());
    var $unitTemperature = $("#unit-temp"); $baseTemperature = parseInt($unitTemperature.val());
  
      $(document).ready(function(e){
        
        $('#id_compound_no').add('#id_kij').change(function(e){
            e.preventDefault()
            var $compNo = $('#id_compound_no').val()
            var $kijMethod = $('#id_kij').val()
            var $compArray = $("#id_chemical0, #id_chemical1, #id_chemical2, #id_chemical3")
            var $fracArray = $("#id_fraction0, #id_fraction1, #id_fraction2, #id_fraction3")
            var $cards = $("#comp1, #comp2, #comp3, #comp4")
            var $bipCards = $("#bip01, #bip02, #bip12, #bip03, #bip13, #bip23")
            var $kijInput = $("#id_kij_input01, #id_kij_input02, #id_kij_input12, #id_kij_input03, #id_kij_input13, #id_kij_input23")
            var $tune = $("#id_tune01, #id_tune02, #id_tune12, #id_tune03, #id_tune13, #id_tune23")
           
            if ($kijMethod == "default"){
              $bipCards.hide()
              $kijInput.attr({"disabled":true, "required":false}); $tune.attr({"disabled":true,"required":false})
              $compArray.attr({"disabled":false,"required":true, "style":"display: block;"}); 
              $fracArray.attr({"disabled":false,"required":true, "style":"display: block;"}); 
              $cards.show();
              $cards.slice(parseInt($compNo), 4).hide();
              $compArray.slice(parseInt($compNo), 4).attr({"disabled":true,"required":false});
              $fracArray.slice(parseInt($compNo), 4).attr({"disabled":true,"required":false});
     
            }
            
            else if($kijMethod == "kij_input"){
              $compArray.attr({"disabled":false,"required":true, "style":"display: block;"});
              $fracArray.attr({"disabled":false,"required":true, "style":"display: block;"});
              $cards.show();
              $cards.slice(parseInt($compNo), 4).hide();
              $compArray.slice(parseInt($compNo), 4).attr({"disabled":true,"required":false});
              $fracArray.slice(parseInt($compNo), 4).attr({"disabled":true,"required":false});
              
              switch(parseInt($compNo)){
                case 1:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  break;
                case 2:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  $kijInput.slice(0,1).attr({"disabled":false,"required":true,"style":"display: block;"});
                  $bipCards.slice(0,1).show()
                  break;
                case 3:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true, "required":false});
                  $kijInput.slice(0,3).attr({"disabled":false,"required":true,"style":"display: block;"});
                  $bipCards.slice(0,3).show()
                  break;
                case 4:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  $kijInput.attr({"disabled":false,"required":true,"style":"display: block;"});
                  $bipCards.show()
                  break;
              }
            }
            else if($kijMethod == "tune"){
              $compArray.attr({"disabled":false,"required":true,"style":"display: block;"});
              $fracArray.attr({"disabled":false,"required":true,"style":"display: block;"});
              $cards.show();
              $cards.slice(parseInt($compNo), 4).hide();
              $compArray.slice(parseInt($compNo), 4).attr({"disabled":true,"required":false});
              $fracArray.slice(parseInt($compNo), 4).attr({"disabled":true,"required":false});
              
              switch(parseInt($compNo)){
                case 1:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  break;
                case 2:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  $tune.slice(0,1).attr({"disabled":false,"required":true,"style":"display: block;"});
                  $bipCards.slice(0,1).show()
                  break;
                case 3:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  $tune.slice(0,3).attr({"disabled":false,"required":true,"style":"display: block;"});
                  $bipCards.slice(0,3).show()
                  break;
                case 4:
                  $bipCards.hide();
                  $kijInput.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":true,"required":false});
                  $tune.attr({"disabled":false,"required":true,"style":"display: block;"});
                  $bipCards.show()
                  break;
              }
            }
        })
  
          var $fugForm = $('.fug-ajax-form');
          $fugForm.submit(function(event){
              event.preventDefault()
              var $myForm = $fugForm.clone(true);
              var entryPressure = $($myForm[0].pressure).val();
              var entryTemperature = $($myForm[0].temperature).val();
              
              $($fugForm[0].pressure).val(convertPressure( parseFloat(entryPressure), parseInt($('#unit-pressure').val()), 0));
              $($fugForm[0].temperature).val(convertTemperature( parseFloat(entryTemperature), parseInt($('#unit-temp').val()), 0));
              var $formData = $fugForm.serialize();
              var $thisURL = $fugForm.attr('data-url') 
              
              $.ajax({
              method: "GET",
              url: $thisURL,
              data: $formData,
              success: function(response){
                $('#result').html(response);
                
                $('html, body').stop().animate({
                  scrollTop: $($('#result-card')).offset().top
              }, 800, 'linear');
              },
          })
          $($fugForm[0].pressure).val(entryPressure)
          $($fugForm[0].temperature).val(entryTemperature)
      })
  
          $unitPressure.change(function(e){
            e.preventDefault();
            var $newUnit = parseInt($("#unit-pressure").val());
            var value = parseFloat($("#id_pressure").val());
            var $newPressure = convertPressure(value, $basePressure, $newUnit).toFixed(3);
            var min = convertPressure(parseFloat($("#id_pressure").attr('min')), $basePressure, $newUnit);
            var max = convertPressure(parseFloat($("#id_pressure").attr('max')), $basePressure, $newUnit);
            $("#id_pressure").attr({"min":min,"max":max});
            $("#id_pressure").prop("value",String($newPressure));
            $basePressure = $newUnit
          })

          $unitTemperature.change(function(e){
            e.preventDefault();
            var $newUnit = parseInt($("#unit-temp").val());
            var value = parseFloat($("#id_temperature").val());
            var $newTemperature = convertTemperature(value, $baseTemperature, $newUnit).toFixed(2);
            var min = convertTemperature(parseFloat($("#id_temperature").attr('min')), $baseTemperature, $newUnit);
            var max = convertTemperature(parseFloat($("#id_temperature").attr('max')), $baseTemperature, $newUnit);
            $("#id_temperature").attr({"min":min,"max":max});
            $("#id_temperature").prop("value",$newTemperature);
            $baseTemperature = $newUnit
          })
  })