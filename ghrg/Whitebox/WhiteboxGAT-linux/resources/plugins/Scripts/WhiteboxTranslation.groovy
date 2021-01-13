import java.text.*
//import java.util.ResourceBundle

class ComparableLocale implements Comparable<ComparableLocale> {
	public String country
	public String language

	public ComparableLocale(Locale l) {
		this.country = l.getDisplayCountry()
		this.language = l.getDisplayLanguage()
	}

	@Override
    public int compareTo(ComparableLocale other) {

    	int ret = this.language.compareTo(other.language)
    	if (ret == 0) {
    		ret = this.country.compareTo(other.country)
    	}
    	return ret
    }
}

println "List of all locales:"

def list = Locale.getAvailableLocales()
def localeList = new ArrayList<ComparableLocale>()
for (Locale aLocale : list) {
    if (aLocale.getDisplayLanguage() != null && !aLocale.getDisplayLanguage().isEmpty()) {
    	localeList.add(new ComparableLocale(aLocale))
    }
}

Collections.sort(localeList)
for (ComparableLocale l : localeList) {
	if (l.country != null && !l.country.isEmpty()) {
		println("${l.language} (${l.country})");
	} else {
		println("${l.language}")
	}
}
        
//def bundle = pluginHost.getMessageBundle()
//println ""
//println bundle.getLocale()
//println ""
//(bundle.getKeys()).each() { it -> println bundle.getObject(it) }
