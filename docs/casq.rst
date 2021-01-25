What CaSQ does
==============

CaSQ [1]_, is a tool that can convert a molecular interaction map built with
CellDesigner_ to an executable Boolean model. The tool is developed in
Python and uses as source the xml file of CellDesigner, in order to infer
preliminary Boolean rules based solely on network topology and semantic
annotations (e.g., certain arcs are noted as catalysis, inhibition, etc.). The
aim is to convert a Process Description representation, i.e., a reaction
model, into a full logical model. The resulting structure is closer to an
Activity Flow diagram, though not in a strict SBGN-PD to SBGN-AF
notion. Moreover logical rules that make the model executable are also
obtained. CaSQ is being used by the `Covid-19 DiseaseMaps consortium`__ to
automatically obtain logical models from maps [2]_.

The conversion happens in 4 steps.

1. First, the map is reduced through a pass of graph-rewriting rules.
   These rules are executed in order and in a single pass, so the rewriting is
   terminating and confluent.

   The idea of this reduction is that a single qualitative species of the
   logical model often represents by its state (active/inactive) several
   species of the original map.

   The rules are the following:

   #. if two species of the map are only reactants in a single reaction, i.e.,
      do not take part in any other reaction, if that reaction is annotated
      as *heterodimer association*, and if one of the reactants is
      annotated as a *receptor*, then the receptor is deleted from the map
      (its annotations are added to the product of the reaction);

   #. if two species of the map take part in a reaction annotated as
      *heterodimer association*, if none of them is annotated as *receptor*,
      and if both do not take active part (i.e., reactant or modifier) in any
      other reaction, then both are merged into the complex, product of the
      reaction (their annotations are added to the product, and the reactions
      that had them as product are rewired to have the complex as product);

   #. if one species only appears in a single reaction as reactant, if that
      reaction has a single product, and if both the reactant and the product
      have the same *name*, then the reactant is deleted (its annotations are
      merged into those of the product);

   #. if one species only appears as reactant in a single reaction (but maybe
      appears as product in another reaction) that has a single product and is
      annotated as *transport*, and if both the reactant and the product have
      the same *name*, then the reactant is merged into the product (its
      annotations are merged into those of the product, and the reactions
      producing it are rewired to the product).

   The rationale of using the *name* to identify the same component in
   different states (gene, RNA, protein, transported/phosphorylated/methylated
   protein, etc.) is that relying on the *active* annotation in CellDesigner
   maps proved to be insufficient and that the model should only keep what
   really contributes to some form of signal propagation.

2. Then, the topology of the model is computed as some simple form of SBGN-PD
   to SBGN-AF conversion, with one qualitative species corresponding to each
   original map species. This species inherits the original map layout, using
   SBML3 Layout package, and MIRIAM annotations (e.g., PubMed IDs as
   *bqbiol:isDescribedBy*). The annotations have been currently arbitrarily
   associated to each regulated component rather than each regulation, but
   that is mostly because tools supporting the latter are quite rare.
   Basically all reactants and modifiers of a reaction get a positive
   influence on all the products of that reaction, whereas all inhibitors get
   a negative influence. Compared to the formal abstraction of influence
   graphs from reaction graphs (Ref FagesSoliman08tcs) note that the mutual
   inhibition between reactants is purposely ignored. This is related to the
   fact that our model simplification step has already condensed active and
   inactive forms of the same species.

3. Finally logical rules of the model are computed. For each species, its
   logical rule is defined as the *OR* over all reactions producing it of
   another *OR* on all positive modifiers (annotated as *catalysis* or
   *transition*) being on, of an *AND* on all products being activated and all
   inhibitors being inactive. Therefore a target is on if one of the reactions
   producing it is on, a reaction is on if all reactants are on, all
   inhibitors are off and one of the catalysts is on.

4. Model cleaning is then done through optional removal of unconnected
   components. From our experience, only keeping the biggest connected
   component is what makes the more sense in a modelling perspective, however
   it is possible to specify a “minimum size” and keep all connected
   components above that size. Names of the qualitative species are also made
   more precise by adding the original *type*/*modifications* of the species
   (e.g., RNA, phosphorylated, etc.) and if there are still homonyms the
   original compartment is added too.

CaSQ generates two output files, the proper logical model encoded in
SBML-qual, a format that is compatible for further analysis with modelling
tools such as GINsim [3]_ or CellCollective [4]_ and a CSV file that contains
information about the names, the logic formulae and the CellDesigner alias.
The second file is mostly for automated treatment.

.. _[1]: https://academic.oup.com/bioinformatics/article/36/16/4473/5836892
.. _CellDesigner: http://celldesigner.org/
__ https://covid.pages.uni.lu/
.. _[2]: https://www.biorxiv.org/content/10.1101/2020.10.26.356014v1
.. _[3]: http://ginsim.org/
.. _[4]: https://cellcollective.org/

CaSQ functions
==============

.. automodule:: casq.celldesigner2qual
   :members:
