/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.distributions;

/**
 *
 * @author Vitaly
 */
public class ConstDistr implements Distribution {

    @Override
    public double getValue() {
        return 1;
    }

    @Override
    public String getName() {
        return "ConstDistr";
    }
    
}
