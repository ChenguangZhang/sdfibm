#pragma once
#include <unordered_map>
#include <iostream>
#include "types.h"
#include "genericfactory.h"

namespace sdfibm
{
    template <typename BaseType, typename ParaType=Foam::dictionary>
    class EntityLibrary : public std::unordered_map<std::string, std::unique_ptr<BaseType>>
    {
    public:
        EntityLibrary() = default;
        EntityLibrary(const dictionary& def)
        {
            for (int i=0; i < def.size(); ++i)
            {
                const dictionary& def_ = def.subDict(def.toc()[i]);
                std::string type = Foam::word(def_.lookup("type"));
                std::string name = Foam::word(def_.lookup("name"));

                this->emplace(name, std::move(GenericFactory<BaseType, ParaType>::create(type, def_)));
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const EntityLibrary<BaseType>& lib)
        {
            os << "Objects in the library:\n";
            int i = 1;
            for (const auto& [key, value] : lib)
            {
                (void) value;
                os << '[' << i++ << "] " << key << std::endl;
            }
            return os;
        }
    };
}