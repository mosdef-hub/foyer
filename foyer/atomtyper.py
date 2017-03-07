from warnings import warn

from oset import oset as OrderedSet

from foyer.exceptions import FoyerError
from foyer.rule import Rule

RULE_NAME_TO_RULE = dict()


def find_atomtypes(atoms, forcefield):
    """Determine atomtypes for all atoms.

    This function is fairly general in that it can function on any list of atom
    objects as long as they have a property `neighbors`, which is a list of
    other atoms that they are bonded to as well as the attributes `whitelist`
    and `blacklist` which are sets - if they are ordered sets it simplifies
    debugging a little bit.

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms whose atomtypes you are looking for.
    forcefield : simtk.openmm.app.Forcefield object
        The forcefield object.
    debug : bool, default=True
        Provides debug information about the logical consistency of the
        atomtyping rules.

    See also
    --------
    _sanitize

    """

    for atom in atoms:
        atom.whitelist = OrderedSet()
        atom.blacklist = OrderedSet()

    _load_rules(forcefield)
    _iterate_rules(atoms, max_iter=10)
    _resolve_atomtypes(atoms)


def _load_rules(forcefield):
    global RULE_NAME_TO_RULE
    RULE_NAME_TO_RULE = dict()
    for rule_name, smarts in forcefield._atomTypeDefinitions.items():
        overrides = forcefield._atomTypeOverrides.get(rule_name)
        RULE_NAME_TO_RULE[rule_name] = Rule(rule_name, forcefield.parser, smarts, overrides=overrides)


def _iterate_rules(atoms, max_iter=10):
    """Iteratively run all the rules until the white- and backlists converge.

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms whose atomtypes you are looking for.
    max_iter : int, optional, default=10
        The maximum number of iterations.

    """
    for _ in range(max_iter):
        max_iter -= 1
        found_something = False
        for atom in atoms:
            for rule in RULE_NAME_TO_RULE.values():
                if rule.name not in atom.whitelist:
                    if rule.matches(atom):
                        atom.whitelist.add(rule.name)
                        atom.blacklist |= rule.overrides
                        found_something = True
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")


def _resolve_atomtypes(atoms):
    """Determine the final atomtypes from the white- and blacklists."""
    for i, atom in enumerate(atoms):
        atomtype = [rule_name for rule_name in atom.whitelist - atom.blacklist]

        if len(atomtype) == 1:
            atom.id = atomtype[0]
        elif len(atomtype) > 1:
            raise FoyerError("Found multiple types for atom {0} ({1}): {2}.".format(i, atom.element.name, atomtype))
        else:
            raise FoyerError("Found no types for atom {0} ({1}).".format(i, atom.element.name))
