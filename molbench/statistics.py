from .comparison import Comparison
from . import logger as log
from collections import defaultdict
import numpy


class Statistics:
    """
    Class for statistical evaluation of a data set.
    """

    available_error_measures = {}

    def __init__(self, data: Comparison) -> None:
        if not isinstance(data, Comparison):
            log.error("Data for statistics evaluation has to be provided as "
                      f"{Comparison}.", "Statistics", "TypeError")
        self._data = data

    @property
    def data(self):
        return self._data

    def compare(self, interest: dict, reference: dict, relative: bool = False,
                relative_damping: float = 0.0) -> dict:
        """
        Computes the signed error for a subset of data as
        interest_value - reference_value.
        If relative is true, it instead computes the relative signed errors as
        (interest_value - reference_value) / abs(reference_value)
        The subset of data for wich to compute the signed error has to be
        provided as two descriptive dictionaries that define which values in
        the dataset should be used as reference and interest values. The keys
        in the dictionary have to match the keys used in the Comparison class
        to separate the data points. Note that all values are first probed to
        be a reference value. Therefore providing an empty reference
        description will result in all values being assigned as reference and
        thus no signed errors will be computed in this case. In a second step,
        reference and interest values are mapped onto each other.
        Thereby, we apply the following rules to the provided descriptions:
        1) If a descriptor ('method', 'basis', ...) is given in both
           dictionaries, the corresponding values are fixed for interest and
           reference -> the inpiut data provides a clear mapping that we don't
           interfere with.
        2) If a descriptor is given in one of the dictionaries, it restricts
           the data in one of the subsets but no further restriction is applied
           to the other subset, e.g., defining the 'method' in the reference
           dictionary restricts the reference data to the given method, but
           compares the reference values with all available methods in the
           interest data subset.
        3) If a descriptor is not given in either of the dictionaries, the
           descriptor of interest and reference values have to be equal, e.g.,
           if 'basis' is not given, only values with equal 'basis' descriptor
           will be compared. An exception to the this rule is the 'data_id'
           descriptor. Since 'data_id' is used to avoid conflicting entries in
           the data set, the 'data_id' descriptor is not forced to be the same
           if no 'data_id' is given in the input.
        Each value is only mapped onto a single reference value.
        """

        identifier = self.identify(interest, reference)
        interest_finder = self.get_interest_values(interest, reference)
        return self._compare(identifier, interest_finder, relative, relative_damping)

    def identify(self, interest: dict, reference: dict) -> callable:
        """Returns a callable to identify whether a value is a reference or
           interest value.
        """
        all_separators = self.data.structure

        def _identify(separators) -> str | None:
            local_dict = {k: v for k, v in zip(all_separators, separators)}
            if reference.items() <= local_dict.items():
                return "reference"
            elif interest.items() <= local_dict.items():
                return "interest"
            else:
                return None
        return _identify

    def get_interest_values(self, interest, reference) -> callable:
        """Returns a callable that identifies the interest values that belong
           to a given reference value.
        """
        common_keys = interest.keys() & reference.keys()
        fixed_interest_separators = {k: interest[k] for k in common_keys}
        all_separators = self.data.structure

        def _get_interest_values(ref_separators: list,
                                 interest_pool: list) -> list:
            separators = []
            for ref_sep, sep in zip(ref_separators, all_separators):
                # key already fixed in the input
                if sep in fixed_interest_separators:
                    separators.append(fixed_interest_separators[sep])
                # data_id special case: if not given twice -> allow all
                # or key given once -> allow all
                elif sep == "data_id" or sep in reference:
                    separators.append(None)
                elif sep in interest:  # defined
                    separators.append(interest[sep])
                else:
                    # key not given -> has to be the same as in reference
                    separators.append(ref_sep)
            # filter out all not compatible values
            interest_values = []
            assigned_interest = []
            for i, (keys, value) in enumerate(interest_pool):
                # check if the keys are compatible with the separators
                compatible = True
                for desired_sep, interest_sep in zip(separators, keys):
                    if desired_sep is not None and desired_sep != interest_sep:
                        compatible = False
                        break
                if compatible:
                    interest_values.append((keys, value))
                    assigned_interest.append(i)
            for i in reversed(assigned_interest):
                del interest_pool[i]
            if len(interest_values) > 1:
                log.warning("Found more than 1 interest value for reference "
                            f"value {ref_separators}.", "Statistics: get_interest_values")
            return interest_values
        return _get_interest_values

    def _compare(self, identify: callable,
                 get_interest_values: callable, 
                 relative: bool = False,
                 relative_damping: float = 0.0) -> dict:
        reference = []
        interest = []
        for keys, value in self.data.walk_values():
            role = identify(keys)
            if role is None:  # neither ref nor interest value
                continue
            elif role[0] == "r":
                reference.append((tuple(keys), value))
            elif role[0] == "i":
                interest.append((tuple(keys), value))
            else:
                log.error(f"Could not assign a role to {keys}.", self,
                          ValueError)

        signed_errors = defaultdict(dict)
        for (ref_keys, ref) in reference:
            interest_values = get_interest_values(ref_keys, interest)
            for interest_keys, values in interest_values:
                se = values - ref
                if relative:
                    se /= abs(ref) + relative_damping
                signed_errors[ref_keys][interest_keys] = se
        
        for ref_keys in signed_errors:
            for interest_keys in signed_errors[ref_keys]:
                if ((abs(signed_errors[ref_keys][interest_keys]) > 1 and not relative)
                     or (abs(signed_errors[ref_keys][interest_keys]) > 10.0 and relative)):
                    log.warning(f"Large Error detected: {signed_errors[ref_keys][interest_keys]}\n"
                                + f"Reference:    {ref_keys}\nInterest:     {interest_keys}\n"
                                + f"Relative Error: {relative}\nDamping: {relative_damping}\n"
                                + "Please check that all calculations involved in this "
                                + "calculation were successful.", "Statistics._compare")
        return signed_errors

    def evaluate(self, signed_errors: dict, *statistical_error_measures,
                 assign: callable = None, proptype: str = None) -> dict:
        """
        Evaluates statistical error measures for the given set of
        signed errors. Statistical error measures can be requested by
        their name, e.g., 'mse' for the mean signed error.
        Optionally, a callable can be provided to determine whether a specific
        signed error should be included in the statistical evaluation.
        By default the signed errors are filtered according to the property
        type ('energy', ...), which can be provided as another optional
        argument.
        """
        statistical_error_measures = set(
            measure.lower() for measure in statistical_error_measures
        )
        if "all" in statistical_error_measures:
            statistical_error_measures.update(
                self.available_error_measures.keys()
            )
            statistical_error_measures.remove("all")

        if assign is None:
            if proptype is None:
                log.error("No assign callable or proptype given.", "Statistics: evaluate")
                return
            assign = self.assign_by_proptye(proptype)

        ret = {}
        for error_measure in statistical_error_measures:
            callback = self.available_error_measures.get(error_measure, None)
            if callback is None:
                log.error("Can not evalute the unknown error measure "
                          f"{error_measure}.", "Statistics", "ValueError")
                continue
            ret[error_measure] = callback(signed_errors, assign)
        return ret

    @staticmethod
    def assign_by_proptye(intproptype: str, refproptype: str = None):
        if refproptype is None:
            refproptype = intproptype

        def assign(refkeys: tuple, interestkeys: tuple) -> bool:
            return (
                refkeys[-2] == refproptype and interestkeys[-2] == intproptype
            )
        return assign


def register_as_error_measure(function):
    """Decorator to register a function as error measure for statistical
       data evaluation.
    """
    Statistics.available_error_measures[function.__name__.lower()] = function
    return function


def _collect_errors(signed_errors: dict, assign: callable) -> list:
    return [value for refkeys, interest in signed_errors.items()
            for interestkeys, value in interest.items()
            if assign(refkeys, interestkeys)]


@register_as_error_measure
def mse(signed_errors: dict, assign: callable):
    """Computes the mean signed error."""
    errors = _collect_errors(signed_errors, assign)
    return numpy.array(errors).mean(axis=0), len(errors)


@register_as_error_measure
def sde(signed_errors: dict, assign: callable):
    """Computes the standard deviation."""
    errors = _collect_errors(signed_errors, assign)
    return numpy.array(errors).std(axis=0), len(errors)


@register_as_error_measure
def mae(signed_errors: dict, assign: callable):
    """Computes the mean absolute error."""
    errors = _collect_errors(signed_errors, assign)
    if len(errors) == 0:
        return numpy.float64(0), 0
    else:
        return (
            numpy.sum(numpy.absolute(e) for e in errors) / len(errors),
            len(errors)
        )


@register_as_error_measure
def min(signed_errors: dict, assign: callable):
    """Computes the minimal signed error."""
    errors = _collect_errors(signed_errors, assign)
    return numpy.array(errors).min(axis=0), len(errors)


@register_as_error_measure
def max(signed_errors: dict, assign: callable):
    """Computes the maximal signed error."""
    errors = _collect_errors(signed_errors, assign)
    return numpy.array(errors).max(axis=0), len(errors)


@register_as_error_measure
def median_se(signed_errors: dict, assign: callable):
    """Computes the median signed error."""
    errors = _collect_errors(signed_errors, assign)
    return numpy.median(numpy.array(errors), axis=0), len(errors)

@register_as_error_measure
def rmsd(signed_errors: dict, assign: callable):
    """Computes the root mean-square deviation."""
    errors = _collect_errors(signed_errors, assign)
    return numpy.sqrt(numpy.mean(numpy.square(errors), axis=0)), len(errors)
