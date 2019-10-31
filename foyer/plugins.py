import foyer
import glob

def collect_plugins(plugin_names=None):
    """
    Detects which plugins are installed

    Arguments
    ---------
    plugin_names: list, default=None
        A list of non-standard plugin names (strings) to additionally
        check for.
    """

    plugin_funcs = [func for func in dir(foyer.forcefields)
                    if "load" in func and "__" not in func]

    if plugin_names:
        plugin_funcs += plugin_names

    plugins = dict()
    for plugin_func in plugin_funcs:
        plugin_loader_func = eval("foyer.forcefields.{}".format(plugin_func))
        # Assumes all plugins are named "load_{plugin_name}"
        plugin_name = "_".join(plugin_func.split("_")[1:])
        # TODO: plugin_version = get_version_info
        plugin_dir = eval("foyer.forcefields.{}.__globals__['__file__']".format(plugin_func))
        # TODO: plugin_xml_path = get_xml_path
        # This assumes that plugin directory tree is consistent.
        # Does not consider versioned FFs.
        plugin_xml_path = glob.glob("{}/xml/*xml".format(plugin_dir))

        plugin = {plugin_name : {"version"       : None,
                                 "load_function" : plugin_loader_func,
                                 "xml_path"      : plugin_xml_path
                                 }}
        plugins.update(plugin)

    return plugins
