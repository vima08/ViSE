/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr;

import altr.entity.Person;
import altr.experiment.Stage;
import altr.experiment.Experiment;
import altr.distributions.NormDistr;
import altr.managers.experiment.AltrExperimentManager;
import altr.managers.experiment.HierarchicalExperimentM;
import altr.managers.experiment.SimpleEgoExperimentManager;
import altr.managers.experiment.SimpleExperimentManager;
import altr.managers.experiment.SimpleThresholdExperimentM;
import altr.managers.experiment.SimpleTwoGroupExperimentM;
import altr.managers.experiment.StrangeAltrExperimentManager;
import javenue.csv.Csv;

/**
 *
 * @author Vitaly
 */
public class Runner {  
    public static Csv.Writer writer = new Csv.Writer("result.csv").delimiter(';');
    
    static int ITERATION_NUMBER = 1000;
    static int STEP_NUMBER = 40;
    static int PEOPLE_NUM = 100;
    
    public static void main(String[] args) throws CloneNotSupportedException {         

       // Person groupMan = new Person(10, true, 0, 1, 0, "name", null);
        Person egoist = new Person(10, true, 1, 0, 0, "name", null);
        Person altrMan = new Person(10, true, 0, 0, 1, "name", null);
                
        //for(int i = 1; i < PEOPLE_NUM - 50; i++){
//        for(int i = 1; i < PEOPLE_NUM; i++){
//        //for(int t = 1; t < 400; t++){
//            Experiment exp = new Experiment(1, ITERATION_NUMBER, "", "test1");
//            Stage stage = new Stage(exp.getId(), STEP_NUMBER, new NormDistr(-.1, 1), 1);
//            exp.addStage(stage);
//            //System.out.println(exp);
//            System.out.println("group size = " + i);
//            //System.out.println("t = " + t);
//            Environment env = new Environment("test");
//            SimpleExperimentManager eM = new SimpleExperimentManager(exp, env, groupMan, 0, egoist, PEOPLE_NUM - 0);
//            //SimpleTwoGroupExperimentM eM = new SimpleTwoGroupExperimentM(exp, env, groupMan, i, egoist, PEOPLE_NUM - i - 50);
//            //AltrExperimentManager eM = new AltrExperimentManager(exp, env, groupMan, i, egoist, PEOPLE_NUM - i);
//            //SimpleThresholdExperimentM eM = new SimpleThresholdExperimentM(exp, env, groupMan, egoist, PEOPLE_NUM, t * 0.001);
//            eM.carryOut();
//        }
  
        /// fdsghdfhd
        
//        for(int i = 1; i < 100; i++){
//            Experiment exp = new Experiment(1, ITERATION_NUMBER, "", "test1");
//            Stage stage = new Stage(exp.getId(), STEP_NUMBER, new NormDistr(-.3, 1), 1);
//            exp.addStage(stage);
//            Environment env = new Environment("test");
//            HierarchicalExperimentM eM = new HierarchicalExperimentM(exp, env, i * 0.01);
//            eM.carryOut();    
//        }
        
//        for(int i = 1; i < 100; i++){
//            Experiment exp = new Experiment(1, ITERATION_NUMBER, "", "test1");
//            Stage stage = new Stage(exp.getId(), STEP_NUMBER, new NormDistr(-5, 10), 1);
//            exp.addStage(stage);            
//            System.out.println("alpha = " + i);            
//            Environment env = new Environment("test");
//            SimpleEgoExperimentManager eM = new SimpleEgoExperimentManager(exp, env, egoist, PEOPLE_NUM, ((double)i )/ 100);
//            eM.carryOut();
//        }
        
//        for(int i = 1; i < 100; i++){
//            Experiment exp = new Experiment(1, ITERATION_NUMBER, "", "test1");
//            Stage stage = new Stage(exp.getId(), STEP_NUMBER, new NormDistr(-.3, 7), 1);
//            exp.addStage(stage);            
//            System.out.println("alpha = " + i);            
//            Environment env = new Environment("test");
//            StrangeAltrExperimentManager eM = new StrangeAltrExperimentManager(exp, env, altrMan, i, egoist, PEOPLE_NUM - i, 0.8, 0.25);
//            eM.carryOut();
//        }
        
        for(int i = 1; i < 100; i++){
            Experiment exp = new Experiment(1, ITERATION_NUMBER, "", "test1");
            Stage stage = new Stage(exp.getId(), STEP_NUMBER, new NormDistr(-.3, 7), 1);
            exp.addStage(stage);            
            System.out.println("delta = " + i);            
            Environment env = new Environment("test");
            StrangeAltrExperimentManager eM = new StrangeAltrExperimentManager(exp, env, altrMan, 30, egoist, PEOPLE_NUM - 30, ((double)i )/ 100, 0.4);
            eM.carryOut();
        }
        
        writer.close();
    }
}
