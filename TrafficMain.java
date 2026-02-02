import java.util.Arrays;

/*
 * COMP2230 - Algorithms
 * Assignment 1 - Main
 * @author  Hakim Newton
 * 
 * This file runs the assignment; you may modify this for testing purposes, but for marking you must ensure that your code works with the original version of this file.
 * The marker will add their own code under the test comment indicated below.
 */

public class TrafficMain {

    int RNG_SEED = 0; // Set this to whatever you like, this ensures that your results are reproducible across runs
    public static void main(String[] args) {
        TrafficMain program = new TrafficMain();
        program.run(args);
    }

    public void run(String[] args){
        TrafficAnalyser analyser = new TrafficAnalyser(RNG_SEED);
        TrafficCheck verifier = new TrafficCheck();
        // You may run the program using the test map by passing "test" as an argument
        if(args.length > 0 && args[0].equals("test")){
            analyser.cityMap = TESTMAP;
        }
        
        analyser.loadMap();
        verifier.storeVerifyMap(analyser.cityMap);

        if(args.length > 0 && args[0].equals("test")){
            checkTestmap(analyser, verifier);
        }

        // Test your calling the methods you have implemented in TrafficAnalyser
        // Use TrafficVerifier to check your results.
        // The following additional methods are available in TrafficCheck to help you with testing:
        // String getRandomIntersection() // Call storeVerifyMap first, will return a random intersection name
        // String getRandomRoad() // Call storeVerifyMap first, will return a random road name
        // Your test code may be placed below this line:

        
    }

    // The following test map is shown in the assignment specification (Figure 1), you may use this to test your code
    // Sample test cases are provided below
    private final String TESTMAP = 
        "{{Red Street, North Crossing, South Junction, 2}, {Yellow Street, South Junction, Western Metro, 9}, "
      + "{Green Street, South Junction, Eastern Plaza, 6}, {Blue Street, Eastern Plaza, Western Metro, 4}, "
      + "{Orange Street, Eastern Plaza, North Crossing, 2}, {Purple Street, Tall Towers, Simple Shops, 4}, "
      + "{Pink Road, Interesting Intersection, Simple Shops, 3}, {Cyan Road, Central Station, Big Boulevard, 7}, "
      + "{Mauve Road, Tall Towers, Perfect Park, 2}, {Violet Road, Perfect Park, Central Station, 8}, "
      + "{Indigo Road, Perfect Park, Big Boulevard, 1}, {Lime Avenue, Perfect Park, Tiny Trainyard, 4}, "
      + "{Gold Avenue, Tiny Trainyard, Big Boulevard, 6}, {Magenta Avenue, Central Station, Tiny Trainyard, 9}}";

    /*
     * Test cases for the test map provided above:
     * 
     * isInInnerCity:
     * North Crossing -> false
     * Eastern Plaza -> false
     * Western Metro -> false
     * Simple Shops -> true
     * Perfect Park -> true
     * Big Boulevard -> true
     * 
     * countInnerCitySlowRoads:
     * Threshold (5) -> 4
     * Threshold (7) -> 2
     * Threshold (4) -> 4
     * 
     * cityBottleneckRoads: (Order is not important)
     * [Pink Road, Purple Street, Mauve Road]
     * 
     * lockdownIntersection:
     * Perfect Park, 1 -> [Mauve Road, Gold Avenue, Purple Street, Lime Avenue, Magenta Avenue, Violet Road, Indigo Road, Cyan Road]
     * Simple Shops, 2 -> [Pink Road, Mauve Road, Purple Street, Lime Avenue, Violet Road, Indigo Road]
     * North Crossing, 4 -> [Orange Street, Green Street, Blue Street, Yellow Street, Red Street]
     * 
     * containProtests:
     * [North Crossing, Simple Shops, Tall Towers] -> 5
     * [South Junction, Eastern Plaza, Interesting Intersection, Central Station, Big Boulevard] -> 5
     * [Tall Towers, Central Station] -> 6
     * 
     */

    // The following method is used to check the test map. It also illustrates the methods available in TrafficVerifier
    private void checkTestmap(TrafficAnalyser analyser, TrafficCheck verifier) {

        System.out.print("CHECK: isInInnerCity: North Crossing (false) ");
        if (analyser.isInInnerCity("North Crossing") ==
            verifier.verify_isInInnerCity("North Crossing"))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: isInInnerCity: Eastern Plaza (false) ");
        if (analyser.isInInnerCity("Eastern Plaza") ==
            verifier.verify_isInInnerCity("Eastern Plaza"))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: isInInnerCity: Western Metro (false) ");
        if (analyser.isInInnerCity("Western Metro") ==
            verifier.verify_isInInnerCity("Western Metro"))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: isInInnerCity: Simple Shops (true) ");
        if (analyser.isInInnerCity("Simple Shops") ==
            verifier.verify_isInInnerCity("Simple Shops"))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: isInInnerCity: Perfect Park (true) ");
        if (analyser.isInInnerCity("Perfect Park") ==
            verifier.verify_isInInnerCity("Perfect Park"))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: isInInnerCity: Big Boulevard (true) ");
        if (analyser.isInInnerCity("Big Boulevard") ==
            verifier.verify_isInInnerCity("Big Boulevard"))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: countInnerCitySlowRoads: Threshold (5) -> (4) ");
        if (analyser.countInnerCitySlowRoads(5) ==
            verifier.verify_countInnerCitySlowRoads(5))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: countInnerCitySlowRoads: Threshold (7) -> (2) ");
        if (analyser.countInnerCitySlowRoads(7) ==
            verifier.verify_countInnerCitySlowRoads(7))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: countInnerCitySlowRoads: Threshold (4) -> (4) ");
        if (analyser.countInnerCitySlowRoads(4) ==
            verifier.verify_countInnerCitySlowRoads(4))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: cityBottleneckRoads: ");
        String [] myAnswer = analyser.cityBottleneckRoads();
        String [] solution = verifier.verify_cityBottleneckRoads();
        // Order is not important, so we sort both arrays before checking them to ensure they will be equal if the same elements are present
        Arrays.sort(myAnswer);
        Arrays.sort(solution);
        if (Arrays.equals(myAnswer,
                          solution))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: lockdownIntersection: Perfect Park, 1 ");
        myAnswer = analyser.lockdownIntersection("Perfect Park", 1);
        solution = verifier.verify_lockdownIntersection("Perfect Park", 1);
        Arrays.sort(myAnswer);
        Arrays.sort(solution);
        if (Arrays.equals(myAnswer,
                          solution))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: lockdownIntersection: Simple Shops, 2 ");
        myAnswer = analyser.lockdownIntersection("Simple Shops", 2);
        solution = verifier.verify_lockdownIntersection("Simple Shops", 2);
        Arrays.sort(myAnswer);
        Arrays.sort(solution);
        if (Arrays.equals(myAnswer,
                          solution))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: lockdownIntersection: North Crossing, 4 ");
        myAnswer = analyser.lockdownIntersection("North Crossing", 4);
        solution = verifier.verify_lockdownIntersection("North Crossing", 4);
        Arrays.sort(myAnswer);
        Arrays.sort(solution);
        if (Arrays.equals(myAnswer,
                          solution))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: containProtests: [North Crossing, Simple Shops, Tall Towers] ");
        if (analyser.containProtests(new String[]{"North Crossing", "Simple Shops", "Tall Towers"}) ==
            verifier.verify_containProtests(new String[]{"North Crossing", "Simple Shops", "Tall Towers"}))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }

        System.out.print("CHECK: containProtests: [Tall Towers, Central Station] ");
        if (analyser.containProtests(new String[]{"Tall Towers", "Central Station"}) ==
            verifier.verify_containProtests(new String[]{"Tall Towers", "Central Station"}))
        {
            System.out.print("-> Passed\n");
        } else {
            System.out.print("-> Failed\n");
        }
    }
}