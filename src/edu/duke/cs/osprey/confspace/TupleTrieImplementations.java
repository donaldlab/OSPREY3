package edu.duke.cs.osprey.confspace;

import java.util.List;

public class TupleTrieImplementations {

    public static class TupETrie extends TupleTrie<TupE> {
        public TupETrie(List<SimpleConfSpace.Position> positions) {
            super(positions);
        }
        @Override
        protected TupE makeT(String repr) {
            return TupE.fromString(repr);
        }
    }

    public static class TupEMappingTrie extends TupleTrie<TupEMapping> {
        public TupEMappingTrie(List<SimpleConfSpace.Position> positions) {
            super(positions);
        }
        @Override
        protected TupEMapping makeT(String repr) {
            return TupEMapping.fromString(repr);
        }
    }

    public static class TupMappingTrie extends TupleTrie<TupMapping>{
        public TupMappingTrie(List<SimpleConfSpace.Position> positions) {
            super(positions);
        }
        @Override
        protected TupMapping makeT(String repr) {
            return TupMapping.fromString(repr);
        }
    }
}
