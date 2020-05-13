package edu.duke.cs.osprey.confspace;

import org.junit.Test;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * Class to test all clases that extend {@link TupE}
 */
public class TestTupE {

    /*
    Testing for TupE base class
     */
    @Test
    public void testTupEFromString(){
        String tupString = "[5=8,6=8,7=4,8=5]->1.4466163083210901";
        TupE test = TupE.fromString(tupString);
        System.out.println(test.toString());
    }

    /*
    Testing for MappableTupE
     */
    @Test
    public void testMappableTupEFromString(){
        String tupString = "[5=8,6=8,7=4,8=5]->[5=124,6=1091,8=0]->1.4466163083210901";
        RCTuple tup = new RCTuple(5,8,6,8,7,4,8,5);
        RCTuple mapTup = new RCTuple(5,124,6,1091,8,0);
        TupEMapping test = TupEMapping.fromString(tupString);

        assertThat(test.toString(), is(tupString));
        assert(test.tup.equals(tup));
        System.out.println(test.mappedTup);
        System.out.println(mapTup);
        assert(test.mappedTup.equals(mapTup));
        assertThat(test.E, is(1.4466163083210901));
    }
}
