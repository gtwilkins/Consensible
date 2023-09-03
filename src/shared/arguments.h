/*
 * Copyright (C) 2017 Glen T. Wilkins <glen.t.wilkins@gmail.com>
 * Written by Glen T. Wilkins
 * 
 * This file is part of the consensible software package <https://github.com/gtwilkins/Consensible>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <vector>
#include <string>

struct Arguments
{
    Arguments( int argc, char** argv );
    bool setBwtPrefix();
    void updateFileIndex();
    std::vector<std::string> inputs_, outputs_, queries_, input_;
    std::string bwtPrefix_, workDir_, outFolder_, outPrefix_;
    std::vector<std::pair<std::string, std::string>> fileIndex_;
    int pInput_;
    bool reindex_, cleanup_, help_, finished_;
private:
    void addInput( std::string fn );
    void checkWorkingDir();
    std::string getCode( int num );
    std::vector<std::string> getFilenameParts( std::string filestr );
    std::vector<std::string> getFilenames( std::string filestr );
    void getFilenames( std::string filestr, std::vector<std::string>& parts, std::vector<std::string>& filenames, int i );
    void setOutFolder( std::string folder );
    void setOutputs();
    void error( std::string msg );
};


#endif /* ARGUMENTS_H */

